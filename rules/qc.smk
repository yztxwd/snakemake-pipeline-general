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
        mv output/qc/fastqc/$(basename {input} .fastq.gz)_fastqc.html {output.html} 2> /dev/null || true
        mv output/qc/fastqc/$(basename {input} .fastq.gz)_fastqc.zip {output.zip} 2> /dev/null || true
        """

rule multiqc:
    input:
        ["output/qc/fastqc/" + str(i).replace('.fq.gz', '').replace('.fastq.gz', '') + "_fastqc.html" for i in list(samples[["fq1", "fq2"]].values.flatten()) if not pd.isnull(i)]
    output:
        report(directory("output/qc/multiqc"), caption="../report/multiqc.rst", htmlindex="multiqc.html", category="QC")
    params:
        config["multiqc"]["params"],
        fastqc_dir="output/qc/fastqc",
    log:
        "logs/multiqc/multiqc.log"
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        mkdir -p {output}
        multiqc {params} --force \
          -o {output} \
          -n "multiqc.html" \
          {params.fastqc_dir} \
          &> {log}
        """
        

rule count_size:
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
