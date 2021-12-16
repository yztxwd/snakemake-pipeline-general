rule fastq_dump:
    output:
        "data/{sra, SRR[0-9]*\.f[ast]*q\.gz}"
    shell:
        """
        fastq-dump --split-3 -O data/ {wildcards.sra}
        """