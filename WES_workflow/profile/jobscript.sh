#!/bin/bash
set -e
# properties = {properties}
eval "$(/gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/software/conda/bin/conda shell.bash hook)"
conda activate gatk
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/software/jdk1.8.0_271.sh
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/software/R-3.6.1.sh
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/software/samtools-1.11.sh 
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/software/bcftools-1.11.sh 
source /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/software/gatk-4.1.9.0.sh
#conda env export

{exec_job}
