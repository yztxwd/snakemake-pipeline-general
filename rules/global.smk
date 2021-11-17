def is_single_end(sample, rep, unit):
    return pd.isnull(samples.loc[(sample, rep, unit), "fq2"])

def checkcontrol(samples):
    return 'control' in samples['condition'].values

snake_dir = ".."    # define relative directory path

rule samtools_sort_name:
    input:
        "{header}.bam"
    output:
        temp("{header}.sortName.bam")
    params:
        "-n"
    threads:
        config['threads']
    resources:
        cpus=config['threads'],
        mem=config['mem']
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        samtools sort {params} -@ {threads} -o {output} {input}
        """

rule samtools_sort_coord:
    input:
        "{header}.bam"
    output:
        "{header}.sort.bam"
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
        samtools sort {params} -@ {threads} -o {output} {input}
        """

rule samtools_index:
    input:
        "{header}.bam"
    output:
        "{header}.bam.bai"
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
        samtools index -@ {threads} {params} {input} {output}
        """
