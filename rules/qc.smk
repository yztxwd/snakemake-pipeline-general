import os

rule multiqc:
    input:
        # how to integrate both trimmomatic and fastp?
        trimmer=expand("output/trimmed/{sample}-{rep}-{unit}.{type}.fastp.json" if config["trimmer"]=="fastp" else "", zip,
                        sample=samples["sample"], rep=samples["rep"], unit=samples["unit"], 
                        type=["se" if pd.isnull(i) else "pe" for i in samples["fq2"]]),
        aligner=expand("logs/bowtie2/{sample}-{rep}-{unit}.{type}.log" if config["aligner"]=="bowtie2" else "", zip, 
                        sample=samples["sample"], rep=samples["rep"], unit=samples["unit"],
                        type=["se" if pd.isnull(i) else "pe" for i in samples["fq2"]]),
        markDuplicates=expand("output/picard/markDuplicates/{sample}-{rep}-{unit}.markDuplicates.txt", zip, 
                        sample=samples["sample"], rep=samples["rep"], unit=samples["unit"]),
        flagstat=expand("output/mapped/{sample}-{rep}-{unit}.flagstat", zip, 
                        sample=samples["sample"], rep=samples["rep"], unit=samples["unit"]),
        fragmentSize=[f"output/qc/CollectInsertSizeMetrics/{row.sample}-{row.rep}.insert_size_metrics.txt" for row in samples.itertuples() if not pd.isnull(row.fq2)],
        fastqc=["output/qc/fastqc/" + os.path.basename(str(i)).replace('.fastq', '').replace('.fq.gz', '').replace('.fastq.gz', '') + "_fastqc.html" for i in list(samples[["fq1", "fq2"]].values.flatten()) if not pd.isnull(i)]
    output:
        html="output/qc/multiqc/multiqc.html",
        dirname=report(directory("output/qc/multiqc"), caption="../report/multiqc.rst", htmlindex="multiqc.html", category="QC")
    params:
        extra=config["multiqc"]["params"],
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "v1.23.4/bio/multiqc"

def find_fastqc_input(wildcards):
    global samples
    fqs = samples[["fq1", "fq2"]].values.flatten()
    fqs = list(fqs[~pd.isnull(fqs)])
    inputs = ["data/" + i for i in fqs if wildcards.sample in i]
    return inputs

rule fastqc:
    input:
        find_fastqc_input
    output:
        html="output/qc/fastqc/{sample}_fastqc.html",
        zip="output/qc/fastqc/{sample}_fastqc.zip"
    params: ""
    log:
        "logs/fastqc/{sample}.fastqc.log"
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        fastqc {params} --quiet \
          --outdir output/qc/fastqc/ {input[0]} \
          > {log}
        """

rule samtools_flagstat:
    input:
        "output/mapped/{sample}-{rep}-{unit}.flag.bam"
    output:
        "output/mapped/{sample}-{rep, [^-]+}-{unit, [^.]+}.flagstat"
    params: ""
    conda:
        f"{snake_dir}/envs/common.yaml"    
    shell:
        "samtools flagstat {input} > {output}"

rule count_size_deeptools:
    input:
        bam="output/mapped/{sample}-{rep}.merge.sort.bam",
        bai="output/mapped/{sample}-{rep}.merge.sort.bam.bai"
    output:
        png=report("output/qc/bamPEFragmentSize/{sample}-{rep, [^-]+}.hist.png", caption="../report/count_size.rst", category="QC")
    params:
        title="{sample}-{rep}",
        extra="--plotFileFormat png"
    log:
        "logs/bamPEFragmentSize/{sample}-{rep}.log"
    threads:
        config["threads"]
    resources:
        cpus=config["threads"],
        mem=config["mem"]
    conda:
        f"{snake_dir}/envs/deeptools.yaml"
    shell:
        "bamPEFragmentSize --bamfiles {input.bam} --histogram {output.png} {params.extra} -T {params.title} -p {threads} > {log}"

rule count_size_picard:
    input:
        bam="output/mapped/{sample}-{rep}.merge.sort.bam",
        bai="output/mapped/{sample}-{rep}.merge.sort.bam.bai"
    output:
        pdf=report("output/qc/CollectInsertSizeMetrics/{sample}-{rep, [^-]+}.insert_size_hist.pdf", caption="../report/count_size.rst", category="QC"),
        txt="output/qc/CollectInsertSizeMetrics/{sample}-{rep, [^-]+}.insert_size_metrics.txt"
    params:
        extra="M=0.05"
    log:
        "logs/bamPEFragmentSize/{sample}-{rep}.log"
    threads:
        config["threads"]
    resources:
        mem=config["mem"]
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        "picard CollectInsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf} {params.extra}"
