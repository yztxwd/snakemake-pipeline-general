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

Workflows supported by this general pipeline:
| name      | Description |
| ----------- | ----------- |
| [chip-seq-standard-pipeline](https://github.com/yztxwd/chip-seq-standard-pipeline) | For ChIP-seq, MNase-seq analysis |
| [atac-seq-standard-pipeline](https://github.com/yztxwd/atac-seq-standard-pipeline) | For ATAC-seq analysis |
| [dna-methylation-pipeline](https://github.com/yztxwd/dna-methylation-pipeline) | For BS-seq, TAB-seq, methylC-seq |
