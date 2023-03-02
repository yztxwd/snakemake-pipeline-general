#!python3

description = """

    Accepts tsv meta info from sra-explorer(https://sra-explorer.info), format it into samples.tsv and downloads.tsv for pipeline
    
    Usage:
		
        $ python tsv_formatter.py meta.tsv

"""

import sys
import pandas as pd

def main():
    meta = pd.read_table(sys.argv[1], sep="\t")
    
    # downloads.tsv
    downloads = meta[["FastQ filename", "FastQ Aspera URL"]]
    downloads.columns = ["fq", "url"]
    downloads.to_csv("downloads.tsv", header=True, index=False, sep="\t")
    
    # samples.tsv
    temp = meta.groupby(["Title", "Accession"]).apply(lambda chunk: [ "aspera/" + i for i in chunk["FastQ filename"]])
    samples = pd.DataFrame(temp.to_list(), index=temp.index).reset_index()
    if samples.shape[1] == 2:
        # single end
        samples.columns = ["sample", "Accession", "fq1"]
        samples["fq2"] = None
    else:
        # pair end
        samples.columns = ["sample", "Accession", "fq1", "fq2"]
    samples['condition'] = "WT"
    samples['rep'] = "rep1"
    samples['unit'] = "tech1"
    samples = samples.loc[:,["sample", "condition", "rep", "unit", "fq1", "fq2"]]
    samples.to_csv("samples.tsv", header=True, index=False, sep="\t")

if __name__ == "__main__":
    main()
