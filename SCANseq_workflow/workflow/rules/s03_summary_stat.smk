import pandas as pd
import yaml

##### load config and sample sheets #####
configfile: "config/config.yaml"
with open(config["libraries"], "r") as f:
    libraries = yaml.load(f, Loader=yaml.Loader)

barcodes = pd.read_csv(config["barcodes"])["barcode"]

with open(config["samples"], "r") as f:
    samples = yaml.load(f, Loader=yaml.Loader)

dir_smk = config["dir_smk"]

rule summary_stat:
    input:
        demultiplexing = "results/01_demultiplexing/{library}/{barcode}.fastq",
        qc = "results/02_qc/{library}_{barcode}/{library}_{barcode}.fastq",
        full_length_cdna = "results/03_full_length/{library}_{barcode}/{library}_{barcode}_full_length.fa",
        full_length_cdna_len = "results/03_full_length/{library}_{barcode}/{library}_{barcode}_full_length_len.txt",
        genome_mapping = "results/04_genome_mapping/{library}_{barcode}/{library}_{barcode}_sorted.bam",
    output:
        summary = "results/summary_stat/{library}_{barcode}/summary_stat.json",
        done = "results/summary_stat/{library}_{barcode}/DONE"
    log: "results/summary_stat/{library}_{barcode}/log"
    threads: 2
    params:
    shell:
        '''
        python workflow/scripts/summary_stat.py --library {wildcards.library} --barcode {wildcards.barcode} --demultiplexing {input.demultiplexing} --qc {input.qc} --full_length_cdna {input.full_length_cdna} --full_length_cdna_len {input.full_length_cdna_len} --genome_mapping {input.genome_mapping} --output {output.summary} 2>{log} &&
        touch {output.done}
        '''


rule summary_stat_merge:
    input:
        expand("results/summary_stat/{library}_{barcode}/summary_stat.json", library = libraries, barcode = barcodes)
    output:
        merge = "results/summary_stat_merge/summary_stat_merge.csv"
    threads: 2
    script:
        "../scripts/summary_stat_merge.py"
