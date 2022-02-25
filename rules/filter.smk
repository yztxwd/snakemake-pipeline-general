def mark_duplicates_find_input(wildcards):
    global samples
    global config
    sx1 = "se" if any(pd.isnull(samples.loc[wildcards.sample, "fq2"])) else "pe"
    sx2 = config["aligner"]
    return f"output/mapped/{{sample}}-{{rep}}-{{unit}}.{sx1}.{sx2}.sort.bam"

rule mark_duplicates:
    input:
        lambda wildcards: mark_duplicates_find_input(wildcards)
    output:
        bam=temp("output/mapped/{sample}-{rep, [^-]+}-{unit, [^.]+}.markDuplicates.bam"),
        metrics=report("output/picard/markDuplicates/{sample}-{rep, [^-]+}-{unit, [^.]+}.markDuplicates.txt", 
            caption="../report/mark_duplicates.rst", category="Filter")
    group: "bam_filter"
    params:
        config["mark_duplicates"]
    log:
        "logs/picard/markDuplicates/{sample}-{rep}-{unit}.log"
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        picard MarkDuplicates \
            I={input}\
            O={output.bam} \
            M={output.metrics} \
            {params} > {log}
        """

rule samtools_view:
    input:
        "output/mapped/{sample}-{rep}-{unit}.markDuplicates.bam"
    output:
        temp("output/mapped/{sample}-{rep, [^-]+}-{unit, [^.]+}.flag.bam")
    group: "bam_filter"
    params:
        lambda wildcards: ((config["samtools_view"]["se"] if is_single_end(**wildcards) 
            else config["samtools_view"]["pe"]))
    threads:
        config['threads']
    resources:
        cpus=config['threads']
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        samtools view {params} -@ {threads} {input} > {output}
        """

rule mapq_filter:
    input:
        "output/mapped/{sample}-{rep}-{unit}.flag.sortName.bam"
    output:
        temp("output/mapped/{sample}-{rep, [^-]+}-{unit, [^.]+}.flag.filtered.bam"),
    group: "bam_filter"
    params:
        lambda wildcards: (config["filter"]["se"] if is_single_end(**wildcards) 
            else config["filter"]["pe"]) 
    conda:
        "../envs/py3.yaml"
    script:
        "../scripts/reads_filter_smk.py"
