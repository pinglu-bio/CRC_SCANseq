import pandas as pd
import yaml

##### load config and sample sheets #####
configfile: "config/config.yaml"
with open(config["libraries"], "r") as f:
    libraries = yaml.load(f, Loader=yaml.Loader)

barcodes = pd.read_csv(config["barcodes"])["barcode"]

dir_smk = config["dir_smk"]

rule s00_raw_data:
    input:
        config["dir_rawdata"] + "/{library}/{library}.fastq.gz"
    output:
        "results/00_raw_data/{library}/{library}.fastq.gz" 
    log:
        "results/00_raw_data/{library}/log" 
    threads: 1
    shell:
        "ln -s {input} {output} 2>{log}"

rule s01_demultiplexing:
    input:
        "results/00_raw_data/{library}/{library}.fastq.gz" 
    output:
        fq = expand("results/01_demultiplexing/{{library}}/{barcode}.fastq", barcode = barcodes),
        done = "results/01_demultiplexing/{library}/DONE"
    log:
        "results/01_demultiplexing/{library}/log" 
    threads: 20
    params:
        barcode_fa = config["barcode_fa"],
        dir_output = config["dir_smk"]+"/results/01_demultiplexing/{library}"
    shell:
        "nanoplexer -t {threads} -b {params.barcode_fa} -s 31 -i -p {params.dir_output} {input} 2>{log} && touch {output.done}"

rule s02_qc:
    input:
        "results/01_demultiplexing/{library}/{barcode}.fastq"
    output:
        fq = "results/02_qc/{library}_{barcode}/{library}_{barcode}.fastq",
        stat = "results/02_qc/{library}_{barcode}/{library}_{barcode}_stat.txt",
        done = "results/02_qc/{library}_{barcode}/DONE"
    log:
        "results/02_qc/{library}_{barcode}/log"
    threads: 8
    params:
        dir_output = "results/02_qc/{library}_{barcode}",
        filename = "{library}_{barcode}"
    shell:
        '''
        NanoFilt --logfile {log} -q 7 -l 100 {input} > {output.fq} &&
        NanoStat --fastq {input} -t {threads} -o {params.dir_output} -n {params.filename}_stat.txt &&
        touch {output.done}
        '''

rule s03_full_length_pychopper: 
# obtain full length cDNA using pychopper; filter read length < 100 and remove gaps using seqkit seq; trim ployA
    input:
        "results/02_qc/{library}_{barcode}/{library}_{barcode}.fastq"
    output:
        temp_full_length_fq = temp("results/03_full_length/{library}_{barcode}/{library}_{barcode}_full_length_temp.fastq"),
        full_length_fa = "results/03_full_length/{library}_{barcode}/{library}_{barcode}_full_length.fa",
        full_length_len = "results/03_full_length/{library}_{barcode}/{library}_{barcode}_full_length_len.txt",
        report_pdf = "results/03_full_length/{library}_{barcode}/pychopper_report.pdf",
        report_tsv = "results/03_full_length/{library}_{barcode}/pychopper_report.tsv",
        done = "results/03_full_length/{library}_{barcode}/DONE"
    log:
        "results/03_full_length/{library}_{barcode}/log"
    threads: 8
    params:
        primers_fa = "config/barcode_96/primers/{barcode}.primers.fa",
        primer_config = "config/barcode_96/primers/primer_config.txt",
        fail = "results/03_full_length/{library}_{barcode}/FAIL"
    shell:
        '''
        pychopper -m edlib -b {params.primers_fa} -c {params.primer_config} -r {output.report_pdf} -S {output.report_tsv} -t {threads} {input} {output.temp_full_length_fq} 2>{log} &&
        seqkit seq -m 100 -g {output.temp_full_length_fq} |awk 'NR%4==1||NR%4==2' - |sed 's/@/>/g' -|trim_isoseq_polyA -i - -t {threads} -G > {output.full_length_fa} 2>>{log} &&
        seqkit fx2tab --length --name {output.full_length_fa} > {output.full_length_len} &&
        touch {output.done} || touch {params.fail}
        '''
