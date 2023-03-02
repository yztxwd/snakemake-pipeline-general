Snakemake general pipeline

@Jianyu Yang, Pennsylvania State University

General pipeline steps including:

* Trimming (trimmomatic or fastp)
* Mapping (bowtie2 or bwa)
* Marking duplicates (picard)
* Filtering samtools flag
* Genome coverage (deeptools)
* MultiQC
* Fragment size histogram (deeptools)