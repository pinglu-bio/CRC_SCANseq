import argparse
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", type=Path, required=True)
    parser.add_argument("-o", "--output", dest="output", type=Path, required=True)
    return parser.parse_args()


def main(args):
    print("Reading input...")
    input = pd.read_table(
        args.input, header=None, dtype={0: str},
        names=["chr", "start", "end", "gene", "score", "strand"]
    )

    print("Summarizing...")
    records = []
    input_gb_gene = {k: v for k, v in input.groupby("gene")}
    for gene, group in tqdm(input_gb_gene.items(), total=len(input_gb_gene)):
        records.append({
            "gene": gene,
            "len": sum(group["end"] - group["start"]),
            "n_exon": group.shape[0],
        })
    summary = pd.DataFrame.from_records(records)
    
    print("Saving output...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output, sep="\t", header=None, index=False)


if __name__ == "__main__":
    main(parse_args())
