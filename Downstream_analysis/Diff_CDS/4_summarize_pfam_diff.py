import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", type=Path, required=True)
    parser.add_argument("-o", "--output", dest="output", type=Path, required=True)
    return parser.parse_args()


def split_pfam(df):
    pattern = re.compile(r'([^\s]+) "([^"]+)";')
    splitted = pd.DataFrame.from_records(np.vectorize(lambda x: {
        key: val for key, val in pattern.findall(x)
    })(df["pfam"]), index=df.index)
    splitted.columns = "pfam_" + splitted.columns
    return df.assign(**splitted)


def main(args):
    print("Reading input...")
    input = pd.read_table(
        args.input, header=None, dtype={0: str}, usecols=[0, 1, 2, 3, 4, 5, 14],
        names=["chr", "start", "end", "gene", "score", "strand", "pfam"]
    )
    input = split_pfam(input)

    print("Summarizing...")
    records = []
    input_gb_gene = {k: v for k, v in input.groupby("gene")}
    for gene, group in tqdm(input_gb_gene.items(), total=len(input_gb_gene)):
        records.append({
            "gene": gene,
            "len": sum(group["end"] - group["start"]),
            "pfam_gene_id": ",".join(sorted(set(group["pfam_gene_id"]))),
            "pfam_transcript_id": ",".join(sorted(set(group["pfam_transcript_id"]))),
        })
    summary = pd.DataFrame.from_records(records)
    summary["n_pfam_gene_id"] = summary["pfam_gene_id"].map(lambda x: len(x.split(",")))
    summary["n_pfam_transcript_id"] = summary["pfam_transcript_id"].map(lambda x: len(x.split(",")))
    
    print("Saving output...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", header=None, index=False)


if __name__ == "__main__":
    main(parse_args())
