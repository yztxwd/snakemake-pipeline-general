#!/usr/bin/python
import sys, pysam, re
import numpy as np

from optparse import OptionParser, IndentedHelpFormatter

class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """reads_filter_smk.py

snakemake version of reads_filter.py
filter reads given the MAPQ threshold or paired-end fragment size

Example:
    rule mapq_filter:
		input:
			"mapped/{sample}-{unit, [^.]+}.F1804.sort.bam"
		output:
			"mapped/{sample}-{unit, [^.]+}.F1804.filtered.bam"
		params:
			"-m pair -f 0 -t 20"
		conda:
			"envs/py3.yaml"
		script:
			"scripts/reads_filter_smk.py"

Dependencies:
    numpy
    pysam

@Copyright 2020, Jianyu Yang, Southern Medical University
"""

parser = OptionParser(description=description, formatter=CustomHelpFormatter())
parser.add_option('-m', '--mode', dest='mode', help="single/pair end mode")
parser.add_option('-f', '--fragment-size', dest='fragment', help="fragment size threshold for paired-end reads, keep reads < threshold")
parser.add_option('-t', '--mapq-threshold', dest='threshold', help="MAPQ threshold to filter reads")
pattern = re.compile("[ ]+")
option, argument = parser.parse_args(pattern.split(str(snakemake.params)))

#Read in snakemake parameters
input = snakemake.input[0]
output = snakemake.output[0]
mode = option.mode
frag_size_threshold = int(option.fragment)
mapq_threshold = int(option.threshold)

#One could look at the file name too see if we're dealing with a SAM or BAM file, but this is simpler
ifile = pysam.AlignmentFile(input, "rb")
ofile = pysam.AlignmentFile(output, "wb", template=ifile)

#Record how many reads filtered during iteration
filter_mapq = 0
filter_size = 0
filter_chrom = 0

#iterate
if mode == 'pair':
	count = 0
	for read in ifile:
		if count % 2 == 0 :
			r1 = read
			count += 1
		elif count % 2 == 1 :
			r2 = read
			count += 1

			# check if read.name is the same, otherwise exit
			if (r1.query_name != r2.query_name):
				raise Exception("Two consecutive reads don't have the same query name! Make sure bam file has been sorted by name!")

			if (r1.mapping_quality < mapq_threshold or r2.mapping_quality < mapq_threshold):
				filter_mapq += 1
				continue
			
			if (r1.reference_name != r2.reference_name):
				filter_chrom += 1
				continue

			r1_coor = r1.reference_start + 1
			r2_coor = r2.reference_start + 1
			left = r1_coor if r1_coor < r2_coor else r2_coor
			right = r2_coor if r1_coor < r2_coor else r1_coor
			read_length = r2.infer_read_length() if r1_coor < r2_coor else r1.infer_read_length()
			if (frag_size_threshold > 0) and ((right - left + read_length) >= frag_size_threshold):
				filter_size += 1
				continue	

			ofile.write(r1)
			ofile.write(r2)

elif mode == 'single':
	for read in ifile:
		if read.mapping_quality < mapq_threshold:
			filter_mapq += 1
			continue
		
		ofile.write(read)
else:
	raise Exception("Please input a valid mode: single/pair")
		
ifile.close()
ofile.close()
print(f"""
Reads filtered due to MAPQ < {mapq_threshold}: {filter_mapq * 2 if mode=="pair" else filter_mapq}
Reads filtered due to fragment size >= {frag_size_threshold}: {filter_size * 2}
Reads filtered due to inconsistent chromosome name: {filter_chrom}
""")


