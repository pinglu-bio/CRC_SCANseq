#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
bam2matrix based on the pileup data
"""

# ---------
# Change Logs:
#
# ---------

__author__ = 'Li Pidong'
__email__ = 'lipidong@126.com'
__version__ = '0.0.1'
__status__ = 'Dev'


import sys
import argparse
import logging
import pdb

import pandas as pd
import pysam
from tqdm import tqdm

# BED/BAM/wig/bedgraph file: 0-based
# VCF/GFF/SAM/Bigwig file: 1-based

def proc_vcf(vcf_file): # read vcf file, each chrom is a key in vcf_dict; each pos is a key in vcf_dict[chrom]; the value corresponding the pos key is a list [ref, alt]
    vcf_dict = {}
    with open(vcf_file) as f:
        for xx in f:
            if xx.startswith('#'):
                continue
            chrom, pos, id_v, ref, alt = xx.rstrip('\n').split('\t')[:5]
            pos = int(pos) -1 # 1-based 转为 0-based
            vcf_dict[(chrom, pos)] = [ref, alt]
    return vcf_dict


def seq_similarity(x, y):
    if len(x) != len(y):
        return 0.0
    count = 0
    for cx, cy in zip(x, y):
        count += cx == cy
    return count / len(x)


def get_snv_info(vcf_dict, samfile, start=None, end=None):
    """

    """
    snv_list = []
    for (chrom, snv_pos), (ref_base, alt_base) in tqdm(vcf_dict.items()):
        if len(ref_base) > 1 and len(alt_base) == 1 and ref_base[0] == alt_base[0]:
            vtype = "del"
            vlen = len(ref_base) - len(alt_base)
        elif len(ref_base) == 1 and len(alt_base) > 1 and ref_base[0] == alt_base[0]:
            vtype = "ins"
            vlen = len(alt_base) - len(ref_base)
        elif len(ref_base) == len(alt_base) and ref_base != alt_base:
            vtype = "snv/dnv"
            vlen = len(alt_base)
        else:
            print("Unrecognized variant type")

        for pileupcolumn in samfile.pileup(contig=chrom, start=snv_pos, end=snv_pos+1,truncate=True):
            pos = pileupcolumn.pos
            ref_count = 0
            alt_count = 0
            other_count = 0

            if vtype == "snv/dnv":
                for pileupread in pileupcolumn.pileups: # 遍历所有reads
                    if pileupread.is_del:
                        continue
                    query_base = pileupread.alignment.query_sequence[
                        pileupread.query_position:(pileupread.query_position + vlen)
                    ]
                    if query_base == ref_base:
                        ref_count += 1
                    elif query_base == alt_base:
                        alt_count += 1 
                    else:
                        other_count += 1
            elif vtype == "del":
                for pileupread in pileupcolumn.pileups: # 遍历所有reads
                    if pileupread.is_del:
                        continue
                    query_base = pileupread.alignment.query_sequence[
                        pileupread.query_position
                    ]
                    if pileupread.indel == 0:  # and query_base == ref_base[0]
                        ref_count += 1
                    elif pileupread.indel == -vlen:  # and query_base == ref_base[0]
                        alt_count += 1
                    else:
                        other_count += 1
            elif vtype == "ins":
                for pileupread in pileupcolumn.pileups: # 遍历所有reads
                    if pileupread.is_del:
                        continue
                    query_base = pileupread.alignment.query_sequence[
                        pileupread.query_position:(pileupread.query_position + pileupread.indel + 1)
                    ]
                    if pileupread.indel == 0:  # and query_base == ref_base
                        ref_count += 1
                    elif pileupread.indel == vlen and seq_similarity(query_base, alt_base) >= 0.7:
                        alt_count += 1
                    else:
                        other_count += 1

            total_count = ref_count + alt_count + other_count
            if total_count != 0:
                snv_list.append({
                    'chrom': chrom, 'pos': snv_pos + 1, 'ref_base': ref_base, 'alt_base': alt_base,
                    'ref_count': ref_count, 'alt_count': alt_count, 'other_count': other_count, 'total_count': total_count
                })
    return snv_list


def get_args():
    parser = argparse.ArgumentParser(prog='bam2matrix', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--vcf', help='vcf', required=True)
    parser.add_argument('--bam', help='bam', required=True)
    parser.add_argument('--output', required=True,
                        help='结果文件')

    parser.add_argument('--sample', help='sample', type=str, required=True)
    parser.add_argument("--quite", help="increase output verbosity",
                        action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()

def main():
    args = get_args()
    bam_file = args.bam
    vcf_file = args.vcf
    output = args.output

    samfile = pysam.AlignmentFile(bam_file, "rb" )
    vcf_dict = proc_vcf(vcf_file)
    snv_list = get_snv_info(vcf_dict, samfile)

    snv_df = pd.DataFrame.from_records(snv_list, columns=[
        "chrom", "pos", "ref_base", "alt_base",
        "ref_count", "alt_count", "other_count", "total_count"
    ])
    del snv_list
    
    snv_df['ref_ratio'] = snv_df['ref_count'] / snv_df['total_count']
    snv_df['alt_ratio'] = snv_df['alt_count'] / snv_df['total_count']
    snv_df['CellName'] = args.sample
    snv_df.to_csv(output, index=False, sep='\t')




if __name__ == '__main__':
    main()
