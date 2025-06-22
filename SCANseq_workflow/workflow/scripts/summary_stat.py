#!/usr/bin/env python
### luping
### 20230107
import os
import gzip
import pysam
import pandas as pd
import argparse
import csv
import json
# import pdb

parser = argparse.ArgumentParser()

summary_dict = {}

parser.add_argument("--library", dest="library", type=str, required=True)
parser.add_argument("--barcode", dest="barcode", type=str, required=True)

parser.add_argument("--demultiplexing", dest="demultiplexing", type=str, required=True)
parser.add_argument("--qc", dest="qc", type=str, required=True)
parser.add_argument("--full_length_cdna", dest="full_length_cdna", type=str, required=True)
parser.add_argument("--full_length_cdna_len", dest="full_length_cdna_len", type=str, required=True)
parser.add_argument("--genome_mapping", dest="genome_mapping", type=str, required=True)

parser.add_argument("--output", dest="output", type=str, required=True)


cmd_args = parser.parse_args()


summary_dict["CellName"] = "_".join([cmd_args.library, cmd_args.barcode])

summary_dict["demultiplexing_reads"] = int(
    os.popen(" ".join(["wc -l",cmd_args.demultiplexing])).read().split()[0]
)//4
summary_dict["qc_reads"] = int(
    os.popen(" ".join(["wc -l",cmd_args.qc])).read().split()[0]
)//4
summary_dict["full_length_cnda_reads"] = int(
    os.popen(" ".join(["wc -l",cmd_args.full_length_cdna])).read().split()[0]
)//2
summary_dict["full_length_cdna_len_mean"] = int(
    os.popen(" ".join(["cut -f 2",cmd_args.full_length_cdna_len,"|awk '{x+=$1; next};END {print int(x/NR)}'"])).read().split()[0]
)
summary_dict["full_length_cdna_len_median"] = int(
    os.popen(" ".join(["cut -f 2",cmd_args.full_length_cdna_len,"|sort |awk '{x[i++]=$1};END {print x[int(i/2)]}'"])).read().split()[0]
)

summary_dict["genome_mapping_primary_reads"] = int(
    os.popen(" ".join(["samtools view -c -F 260",cmd_args.genome_mapping])).read().split()[0]
)

with open(cmd_args.output, "w") as f:
    json.dump(summary_dict, f)

#mydf = pd.DataFrame.from_records(summary_dict)
#mydf.to_csv(cmd_args.output, index=False)


