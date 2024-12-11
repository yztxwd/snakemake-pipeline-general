rule fastp_pe:
    input:
        r1=lambda wildcards: "data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq1"],
        r2=lambda wildcards: "data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq2"],
    output:
        r1=temp("output/trimmed/{sample}-{rep, [^-]+}-{unit}.fastp.1.fq.gz"),
        r2=temp("output/trimmed/{sample}-{rep, [^-]+}-{unit}.fastp.2.fq.gz"),
        r1_unpaired=temp("output/trimmed/{sample}-{rep, [^-]+}-{unit}.fastp.1.unpaired.fq.gz"),
        r2_unpaired=temp("output/trimmed/{sample}-{rep, [^-]+}-{unit}.fastp.2.unpaired.fq.gz"),
        html=report("output/trimmed/{sample}-{rep, [^-]+}-{unit}.pe.fastp.html", caption="../report/fastp.rst", category="fastp"),
        json="output/trimmed/{sample}-{rep, [^-]+}-{unit}.pe.fastp.json"
    log:
        "logs/fastp/{sample}-{rep}-{unit}.fastp.log"
    params:
        extra=config["fastp"]["extra"]
    threads:
        config["threads"]
    resources:
        cpus=config['threads'],
        mem=config['mem']
    conda:
        f"{snake_dir}/envs/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --unpaired1 {output.r1_unpaired} --unpaired2 {output.r2_unpaired} \
          --html {output.html} \
          --json {output.json} \
            {params.extra} &> {log}
        """

rule fastp_se:
    input:
        lambda wildcards: "data/" + samples.loc[(wildcards.sample, wildcards.rep, wildcards.unit), "fq1"]
    output:
        temp("output/trimmed/{sample}-{rep, [^-]+}-{unit, [^.]+}.fastp.fq.gz"),
        html=report("output/trimmed/{sample}-{rep, [^-]+}-{unit}.se.fastp.html", caption="../report/fastp.rst", category="fastp"),
        json="output/trimmed/{sample}-{rep, [^-]+}-{unit}.se.fastp.json"
    log:
        "logs/fastp/{sample}-{rep}-{unit}.fastp.log"
    params:
        extra=config["fastp"]["extra"]
    threads:
        config["threads"]
    resources:
        cpus=config['threads'],
        mem=config['mem']
    conda:
        f"{snake_dir}/envs/fastp.yaml"
    shell:
        """
        fastp -i {input} \
          -o {output} \
          --html {output.html} \
          --json {output.json} \
          {params.extra} &> {log}
        """

