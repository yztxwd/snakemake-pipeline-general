def align_pe_find_input(wildcards):
    global samples
    global config

    trimmer=config["trimmer"]
    if config["skiptrim"]:
        return ["data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq1"],
                "data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq2"]]
    else:
        return [f"output/trimmed/{{sample}}-{{rep}}-{{unit}}.{trimmer}.1.fq.gz",
                f"output/trimmed/{{sample}}-{{rep}}-{{unit}}.{trimmer}.2.fq.gz"]

def align_se_find_input(wildcards):
    global samples
    global config

    trimmer=config["trimmer"]
    if config["skiptrim"]:
        return "data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq1"]
    else:
        return f"output/trimmed/{{sample}}-{{rep}}-{{unit}}.{trimmer}.fq.gz" 
    
rule bowtie2_mapping_pe:
    input:
        lambda wildcards: align_pe_find_input(wildcards)
    output:
        temp("output/mapped/{sample}-{rep}-{unit, [^.]+}.pe.bowtie2.bam")
    log:
        "logs/bowtie2/{sample}-{rep, [^-]+}-{unit}.log"
    params:
        index=lambda wildcards: config["bowtie2"]["index"],
        extra=config["bowtie2"]["extra"]
    threads: 
        config["threads"]
    resources:
        cpus=config["threads"],
        mem=20
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        bowtie2 --threads {threads} {params.extra} \
            -x {params.index} -1 {input[0]} -2 {input[1]} \
            | samtools view -Sbh -o {output} &> {log}
        """

rule bowtie2_mapping_se:
    input:
        lambda wildcards: align_se_find_input(wildcards)
    output:
        temp("output/mapped/{sample}-{rep}-{unit, [^.]+}.se.bowtie2.bam")
    log:
        "logs/bowtie2/{sample}-{rep, [^-]+}-{unit}.log"
    params:
        index=lambda wildcards: config["bowtie2"]["index"],
        extra=config["bowtie2"]["extra"]
    threads: 
        config["threads"]
    resources:
        cpus=config["threads"],
        mem=20
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        bowtie2 --threads {threads} {params.extra} \
            -x {params.index} -U {input} \
            | samtools view -Sbh -o {output} &> {log}
        """  
