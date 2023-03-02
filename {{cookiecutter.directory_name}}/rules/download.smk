rule fastq_dump_single:
    output:
        temp("data/{sra, SRR[0-9]*}.fastq")
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        fasterq-dump -O data/ {wildcards.sra}
        """

rule fastq_dump_pair:
    output:
        r1=temp("data/{sra, SRR[0-9]*}_1.fastq"),
        r2=temp("data/{sra, SRR[0-9]*}_2.fastq")
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        fasterq-dump --split-3 -O data/ {wildcards.sra}
        """

rule aspera:
    output:
        temp("data/aspera/{accession, [SE]RR*[0-9]*[_12]*}.fastq.gz")
    params:
        url=lambda wildcards: downloads.loc[(wildcards.accession + ".fastq.gz"), "url"],
        private_key=config["aspera_private_key"]
    conda:
        f"{snake_dir}/envs/aspera.yaml"
    shell:
        """
        ascp -QT -l 300m -P33001 -i {params.private_key} {params.url} data/aspera/
        """
