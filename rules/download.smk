rule fastq_dump_single:
    output:
        "data/{sra, SRR[0-9]*}.fastq.gz"
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        fastq-dump -O data/ {wildcards.sra}
        gzip data/{wildcards.sra}.fastq
        """

rule fastq_dump_pair:
    output:
        r1="data/{sra, SRR[0-9]*}_1.fastq.gz",
        r2="data/{sra, SRR[0-9]*}_2.fastq.gz"
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        fastq-dump --split-3 -O data/ {wildcards.sra}
        gzip data/{wildcards.sra}_1.fastq
        gzip data/{wildcards.sra}_2.fastq
        """