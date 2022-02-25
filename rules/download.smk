rule fastq_dump_single:
    output:
        temp("data/{sra, SRR[0-9]*}.fastq.gz")
    group: "download"
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        fastq-dump -O data/ {wildcards.sra}
        gzip data/{wildcards.sra}.fastq
        """

rule fastq_dump_pair:
    output:
        r1=temp("data/{sra, SRR[0-9]*}_1.fastq.gz"),
        r2=temp("data/{sra, SRR[0-9]*}_2.fastq.gz")
    group: "download"
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        fastq-dump --split-3 -O data/ {wildcards.sra}
        gzip data/{wildcards.sra}_1.fastq
        gzip data/{wildcards.sra}_2.fastq
        """

rule aspera:
    output:
        temp("data/aspera/{accession, [SE]RR*[0-9]*[_12]*}.fastq.gz")
    group: "download"
    params:
        url=lambda wildcards: downloads.loc[(wildcards.accession + ".fastq.gz"), "url"],
        private_key=config["aspera_private_key"]
    conda:
        f"{snake_dir}/envs/aspera.yaml"
    shell:
        """
        ascp -QT -l 300m -P33001 -i {params.private_key} {params.url} data/aspera/
        """
