import argparse
import os
import re
import shutil
import tempfile
from multiprocessing import Pool
from pathlib import Path
from subprocess import PIPE, Popen

import numpy as np
import pandas as pd
from tqdm import tqdm


BED = ["chr", "start", "end", "name", "score", "strand"]
gtf_gb_tx, pairing_gb_gene, empty_gtf = None, None, None


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pairing", dest="pairing", type=Path, required=True)
    parser.add_argument("-g", "--gtf", dest="gtf", type=Path, required=True)
    parser.add_argument("-o", "--output", dest="output", type=Path, required=True)
    parser.add_argument("-j", "--n-jobs", dest="n_jobs", type=int, default=8)
    return parser.parse_args()


def split_attribute(gtf):
    pattern = re.compile(r'([^\s]+) "([^"]+)";') # 编译
    splitted = pd.DataFrame.from_records(np.vectorize(lambda x: {
        key: val for key, val in pattern.findall(x)
    })(gtf["attr"]), index=gtf.index)
    del gtf["attr"]
    return gtf.assign(**splitted)


def gtf2bed(gtf, name):
    usecols = ["chr", "start", "end", name, "score", "strand"]
    bed = gtf.loc[:, usecols]
    bed.columns = BED
    bed["start"] -= 1
    return bed


def write_bed(bed, path):
    bed.to_csv(path, sep="\t", header=False, index=False)


def read_bed(path):
    return pd.read_table(path, header=None, names=BED)


def sort(bed):
    return bed.sort_values(["chr", "start"])


def diff(a, b):
    cmd1 = ["bedtools", "subtract", "-s", "-a", a, "-b", b]
    cmd2 = ["bedtools", "subtract", "-s", "-a", b, "-b", a]
    with Popen(cmd1, stdout=PIPE) as p1, Popen(cmd2, stdout=PIPE) as p2:
        return pd.concat([
            read_bed(p1.stdout),
            read_bed(p2.stdout),
        ], ignore_index=True)


def merge(i):
    cmd = ["bedtools", "merge", "-s", "-c", "4,5,6", "-o", "distinct", "-i", i]
    with Popen(cmd, stdout=PIPE) as p:
        return read_bed(p.stdout)


def worker(gene):
    try:
        global pairing_gb_gene, gtf_gb_tx, empty_gtf
        tmp = Path(tempfile.mkdtemp(prefix="PAIREDDIFF"))
        gb_cluster = pairing_gb_gene[gene].groupby("cluster")
        try:
            cc_tx = gb_cluster.get_group("Stem/TA-like")["TranName"]
            ne_tx = gb_cluster.get_group("Stem/TA")["TranName"]
        except KeyError:
            return pd.DataFrame(columns=BED)

        cc_gtf = pd.concat([gtf_gb_tx.get(item, empty_gtf) for item in cc_tx], ignore_index=True)
        ne_gtf = pd.concat([gtf_gb_tx.get(item, empty_gtf) for item in ne_tx], ignore_index=True)

        cc_ori = sort(gtf2bed(cc_gtf, name="transcript_name"))
        write_bed(cc_ori, tmp / "cc_ori.bed")
        cc_merge = sort(merge(tmp / "cc_ori.bed")) if cc_ori.shape[0] else cc_ori
        write_bed(cc_merge, tmp / "cc_merge.bed")

        ne_ori = sort(gtf2bed(ne_gtf, name="transcript_name"))
        write_bed(ne_ori, tmp / "ne_ori.bed")
        ne_merge = sort(merge(tmp / "ne_ori.bed")) if ne_ori.shape[0] else ne_ori
        write_bed(ne_merge, tmp / "ne_merge.bed")

        region_diff = diff(tmp / "cc_merge.bed", tmp / "ne_merge.bed")
        region_diff["name"] = gene
    finally:
        shutil.rmtree(tmp)
    return region_diff


def main(args):
    global pairing_gb_gene, gtf_gb_tx, empty_gtf

    print("Reading input...")
    pairing = pd.read_table(args.pairing)
    pairing_gb_gene = {k: v for k, v in pairing.groupby("GeneName")}
    gtf = pd.read_table(
        args.gtf, header=None,
        names=["chr", "start", "end", "score", "strand", "attr"],
        usecols=[0, 3, 4, 5, 6, 8], dtype={0: str}
    )
    gtf = split_attribute(gtf)
    gtf_gb_tx = {k: v for k, v in gtf.groupby("transcript_name")}
    empty_gtf = gtf.iloc[[]]

    print("Computing difference...")
    try:
        region_diff_list = []
        task = pairing_gb_gene.keys()
        with Pool(args.n_jobs) as pool:
            for region_diff in tqdm(pool.imap(worker, task), total=len(task)):
                region_diff_list.append(region_diff)
    except KeyboardInterrupt:
        print("Interrupted by user.")
    finally:
        print("Saving result...")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        region_diff_all = pd.concat(region_diff_list, ignore_index=True)
        write_bed(region_diff_all, args.output)


if __name__ == "__main__":
    main(parse_args())
