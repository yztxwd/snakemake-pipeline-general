rule merge_bam:
    input:
        lambda wildcards: expand("output/mapped/{sample}-{rep}-{unit}.flag.bam" if config['filter']['skip'] else "output/mapped/{sample}-{rep}-{unit}.flag.filtered.bam",
            **wildcards, unit=samples.loc[(wildcards.sample, wildcards.rep),  'unit'])
    output:
        temp("output/mapped/{sample}-{rep, [^.]+}.merge.bam")
    params:
        ""
    threads:
        config['threads']
    resources:
        cpus=config['threads'],
        mem=config['mem']
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        samtools merge -@ {threads} {params} \
            {output} {input}
        """

