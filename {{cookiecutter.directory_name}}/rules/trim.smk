rule trim_pe:
    input:
        r1=lambda wildcards: "data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq1"],
        r2=lambda wildcards: "data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq2"],
    output:
        r1=temp("output/trimmed/{sample}-{rep, [^-]+}-{unit}.trimmomatic.1.fq.gz"),
        r2=temp("output/trimmed/{sample}-{rep, [^-]+}-{unit}.trimmomatic.2.fq.gz"),
        r1_unpaired=temp("output/trimmed/{sample}-{rep, [^-]+}-{unit}.trimmomatic.1.unpaired.fq.gz"),
        r2_unpaired=temp("output/trimmed/{sample}-{rep, [^-]+}-{unit}.trimmomatic.2.unpaired.fq.gz")
    log:
        "logs/trimmomatic/{sample}-{rep}-{unit}.trimmomatic.log"
    params:
        trimmer=["ILLUMINACLIP:" + config["trimmomatic"]["pe_adapter"] + config["trimmomatic"]["adapter_trimmer"], config["trimmomatic"]["trimmer"]],
        extra=""
    threads:
        config["threads"]
    resources:
        cpus=config['threads'],
        mem=config['mem']
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {params.extra} \
          {input.r1} {input.r2} \
          {output.r1} {output.r1_unpaired} \
          {output.r2} {output.r2_unpaired} \
          {params.trimmer} \
          &> {log}
        """

rule trim_se:
    input:
        lambda wildcards: "data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq1"]
    output:
        temp("output/trimmed/{sample}-{rep, [^-]+}-{unit, [^.]+}.trimmomatic.fq.gz")
    log:
        "logs/trimmomatic/{sample}-{rep}-{unit}.trimmomatic.log"
    params:
        trimmer=["ILLUMINACLIP:" + config["trimmomatic"]["se_adapter"] + config["trimmomatic"]["adapter_trimmer"], config["trimmomatic"]["trimmer"]],
        extra=""
    threads:
        config["threads"]
    resources:
        cpus=config['threads'],
        mem=config['mem']
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        trimmomatic SE -threads {threads} {params.extra} \
          {input} {output} \
          {params.trimmer} \
          &> {log}
        """

