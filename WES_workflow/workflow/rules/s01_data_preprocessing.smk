import pandas as pd
##### load config and sample sheets #####
configfile: "config/config.yaml"

sample_table = pd.read_csv(config["sample_table"])
samples = [s for p in pd.concat([sample_table["normal"], sample_table["tumor"]]).str.split(";") for s in p]
patients = [p for p in sample_table["patient"]]

rule s01_00_rawdata:
    input:
        r1 = config["dir_rawdata"] + "{sample}/{sample}_R1.fq.gz",
        r2 = config["dir_rawdata"] + "{sample}/{sample}_R2.fq.gz"
    output:
        r1 = "results/01_00_rawdata/{sample}/{sample}_R1.fq.gz",
        r2 = "results/01_00_rawdata/{sample}/{sample}_R2.fq.gz"
    log:
        "results/01_00_rawdata/{sample}/log"
    threads: 1
    shell:
        '''
        ln -s {input.r1} {output.r1} 2>{log} && 
        ln -s {input.r2} {output.r2} 2>>{log}
        '''

rule s01_01_qc:
    input:
        r1 = rules.s01_00_rawdata.output.r1,
        r2 = rules.s01_00_rawdata.output.r2
    output:
        r1 = "results/01_01_qc/{sample}/{sample}_R1.fq.gz",
        r2 = "results/01_01_qc/{sample}/{sample}_R2.fq.gz",
        html = "results/01_01_qc/{sample}/{sample}_fastp.html",
        json = "results/01_01_qc/{sample}/{sample}_fastp.json"
    log:
        "results/01_01_qc/{sample}/log"
    threads: 4
    shell:
        '''
        fastp -q 20 -u 40 -w {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json}
        '''

rule s01_02_mapping:
    input:
        r1 = rules.s01_01_qc.output.r1,
        r2 = rules.s01_01_qc.output.r2
    output:
        bam = "results/01_02_mapping/{sample}/{sample}_sort.bam"
    log:
        "results/01_02_mapping/{sample}/log"
    threads: 4
    params:
        rg = r"@RG\tID:{sample}\tLB:{sample}\tPL:ILLUMINA\tPM:HISEQ\tSM:{sample}",
        ref = config["file_ref_dna_fa"]
    shell:
        '''
        bwa mem -Y -v 3 -t {threads} -R "{params.rg}" {params.ref} {input.r1} {input.r2} 2>{log} | samtools view -Sb - | samtools sort -o {output} -
        '''

rule s01_03_duprm:
    input:
        bam = rules.s01_02_mapping.output.bam
    output:
        bam = "results/01_03_duprm/{sample}/{sample}_duprm.bam",
        dup_m = "results/01_03_duprm/{sample}/{sample}_duprm_metrics.txt"
    log:
        "results/01_03_duprm/{sample}/log"
    threads: 20
    shell:
        '''
        gatk MarkDuplicatesSpark -I {input.bam} -O {output.bam} -M {output.dup_m} --optical-duplicate-pixel-distance 2500 --create-output-bam-index true --conf 'spark.executor.cores=4' >{log} 2>&1
        '''

rule s01_04_BaseRecalibrator:
    input:
        bam = rules.s01_03_duprm.output.bam
    output:
        txt = "results/01_04_BaseRecalibrator/{sample}/{sample}_BaseRecal.txt"
    log:
        "results/01_04_BaseRecalibrator/{sample}/log"
    threads: 2
    params:
        ref = config["file_ref_dna_fa"],
        dbsnp_vcf = config["dbsnp_vcf"],
        mills_vcf = config["mills_vcf"],
        g1000_vcf = config["g1000_vcf"],
        exon_target = config["exon_target"]
    shell:
        '''
        gatk BaseRecalibrator -I {input.bam} -R {params.ref} --known-sites {params.dbsnp_vcf} --known-sites {params.mills_vcf} --known-sites {params.g1000_vcf} -L {params.exon_target} -O {output.txt} --use-original-qualities 2>{log}
        '''

rule s01_05_ApplyBQSR:
    input:
        bam = rules.s01_03_duprm.output.bam,
        txt = rules.s01_04_BaseRecalibrator.output.txt
    output:
        bam = "results/01_05_ApplyBQSR/{sample}/{sample}_recal.bam"
    log:
        "results/01_05_ApplyBQSR/{sample}/log"
    threads: 2
    params:
        ref = config["file_ref_dna_fa"],
        exon_target = config["exon_target"]
    shell:
        '''
        gatk ApplyBQSR -I {input.bam} -R {params.ref} -L {params.exon_target} -O {output.bam} -bqsr {input.txt} --add-output-sam-program-record --use-original-qualities 2>{log}
        '''
