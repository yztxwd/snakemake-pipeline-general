rule bwa_mapping_pe:
    input:
        lambda wildcards: align_pe_find_input(wildcards)
    output:
        temp("output/mapped/{sample}-{rep}-{unit, [^.]+}.pe.bwa.bam")
    log:
        "logs/bwa/{sample}-{rep, [^-]+}-{unit}.log"
    params:
        index=lambda wildcards: config["bwa"]["index"],
        extra=config["bwa"]["extra"]
    threads: 
        config["threads"]
    resources:
        cpus=config["threads"],
        mem=config['mem']
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        bwa mem -t {threads} {params.extra} \
            {params.index} {input[0]} {input[1]} \
            | samtools view -Sbh -o {output} > {log}
        """

rule bwa_mapping_se:
    input:
        lambda wildcards: align_se_find_input(wildcards)
    output:
        temp("output/mapped/{sample}-{rep}-{unit, [^.]+}.se.bwa.bam")
    log:
        "logs/bwa/{sample}-{rep, [^-]+}-{unit}.log"
    params:
        index=lambda wildcards: config["bwa"]["index"],
        extra=config["bwa"]["extra"]
    threads: 
        config["threads"]
    resources:
        cpus=config["threads"],
        mem=config['mem']
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        bwa mem -t {threads} {params.extra} \
            {params.index} {input} \
            | samtools view -Sbh -o {output} > {log}
        """  
