#!/data1/yztxwd/miniconda3/envs/py3/bin/python
#coding=utf-8

import sys, os, re
import pandas as pd
import matplotlib.pyplot as plt

from collections import defaultdict
from optparse import OptionParser, IndentedHelpFormatter

class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """count_size_smk.py

snakemake version of count_size.py
Count the fragment size frequency given the bam file

Example:
    rule count_size:
        input:
            "mapped/{sample}.merged.coverage.1b.bg"
        output:
            freq="summary/{sample}.size.freq"
            png="summary/{sample}.size.png"
        conda:
            "envs/py3.yaml"
        script:
            "scripts/count_size_smk.py"

Dependencies:
    numpy
    pandas

@Copyright 2020, Jianyu Yang, Southern Medical University
"""

parser = OptionParser(description=description, formatter=CustomHelpFormatter())

# load input and output from snakemake object
coverage = snakemake.input[0]
freq = snakemake.output.freq
png = snakemake.output.png

# count the fragment size frequency
store = defaultdict(lambda: 0)
with open(coverage, 'r') as f:
    for line in f:
        chr, start, end, depth, id = line.split('\t')
        start = int(start); end = int(end)
        store[end-start+1] += 1

# save the fragment size frequency
df = pd.DataFrame(store.values(), store.keys()).sort_index()
df.columns = [coverage]
df.fillna(0)
df.to_csv(freq, sep='\t')

# plot the frequency plot 
plt.figure(figsize=(15, 15))
plt.plot(list(df.index), list(df[coverage]))
plt.savefig(png, dpi=300)
