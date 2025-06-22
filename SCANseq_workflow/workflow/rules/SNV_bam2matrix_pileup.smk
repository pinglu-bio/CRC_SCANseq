import pandas as pd
import yaml

##### load config and sample sheets #####
configfile: "config/config.yaml"

with open(config["samples"], "r") as f:
    samples = yaml.load(f, Loader=yaml.Loader)
dir_smk = config["dir_smk"]
#dir_script = config["dir_scripts"]

rule snv_pileup:
    input:
        "results/04_genome_mapping/{sample}/{sample}_sorted.bam"
    output:
        tsv = "results/SNV_bam2matrix_pileup/{sample}/{sample}_snv.tsv.gz",
        done = "results/SNV_bam2matrix_pileup/{sample}/DONE_pileup"
    log:
        "results/SNV_bam2matrix_pileup/{sample}/log_pileup"
    threads: 1
    params:
        vcf = "/gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/project/ont_iso_crc/2_RNA_CRC/01_RNA_analysis/08_WES/01_WES/01_WES_filter_nochr.vcf",
        exe_bam2matrix = "workflow/scripts/bam2matrix_pileup_based_on_snv.py"
    shell:
        '''
        python {params.exe_bam2matrix} --vcf {params.vcf} --bam {input} --sample {wildcards.sample} --output {output.tsv} 2>{log} &&
        touch {output.done}
        '''
rule rna_snv_merge:
    input:
        expand("results/SNV_bam2matrix_pileup/{sample}/{sample}_snv.tsv.gz", sample = samples)
    output:
        merge = "results/SNV_bam2matrix_pileup_merge/RNA_SNV_merged_1-based.csv"
    threads: 2
    script:
        "../scripts/snv_merge.py"
