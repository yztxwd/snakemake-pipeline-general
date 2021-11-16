# snakemake-pipeline-general

@Author, Jianyu Yang, Pennsylvania State University

General rules could be used for specific purpose

The following files could be expected from this general pipeline:
- output/mapped/{sample}-{rep}.merge.bam
- output/mapped/{sample}-{rep}-{unit}.flag.bam
- output/qc/multiqc/multiqc.html
- output/qc/bamPEFragmentSize/{sample}-{rep}.hist.png
- output/coverage/{sample}-{rep}.bgToBw.bw
- output/coverage/{sample}-{rep}.bamCov.bw
- output/coverage/{sample}-{rep}.bamCompare.bw