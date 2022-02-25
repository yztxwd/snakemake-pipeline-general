# snakemake-pipeline-general

@Jianyu Yang, Pennsylvania State University

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
| Name      | Description |
| ----------- | ----------- |
| [chip-seq-standard-pipeline](https://github.com/yztxwd/chip-seq-standard-pipeline) | For ChIP-seq, MNase-seq analysis |
| [atac-seq-standard-pipeline](https://github.com/yztxwd/atac-seq-standard-pipeline) | For ATAC-seq analysis |
| [dna-methylation-pipeline](https://github.com/yztxwd/dna-methylation-pipeline) | For BS-seq, TAB-seq, methylC-seq |

Change log:
- Feb 25, 2022: Job grouping for cluster execution
- Dec 21, 2021: Signal coverage profile around TSS integrated
- Nov 23, 2021: Snakemake report feature integrated 
- Nov 20, 2021: All snakemake pipelines have been integrated with cookiecutter, to motivates easily deployment of pipeline
