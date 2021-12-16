import pandas as pd
from snakemake.utils import validate, min_version

#### Set minimum snakemake version ####
min_version("5.1.2")

#### Load config and sample sheets ####

configfile: "config.yaml"
#validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index(["sample", "rep", "unit"], drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

#### target rules ####

rule all:
    input:
        "output/qc/multiqc/",
        expand("output/qc/bamPEFragmentSize/{samples}-{rep}.hist.png", zip, samples=samples["sample"], rep=samples["rep"]),
        expand("output/coverage/{samples}-{rep}.bgToBw.bw", zip, samples=samples["sample"], rep=samples["rep"]),
        expand("output/coverage/{samples}-{rep}.bamCov.bw", zip, samples=samples["sample"], rep=samples["rep"]),
        expand('output/coverage/{samples}-{rep}.bamCompare.bw', zip, 
                 samples=samples.loc[samples["condition"]!="control", "sample"] if "control" in samples["condition"].values else [],
                 rep=samples.loc[samples["condition"]!="control", "rep"] if "control" in samples["condition"].values else []),

#### setup singularity ####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

#### setup report ####
report: "report/workflow.rst"

#### load rules ####
include: "rules/global.smk"
include: "rules/download.smk"
include: "rules/qc.smk"
include: "rules/trim.smk"
include: "rules/fastp.smk"
include: "rules/bowtie2.smk"
include: "rules/bwa.smk"
include: "rules/filter.smk"
include: "rules/merge.smk"
include: "rules/coverage.smk"
