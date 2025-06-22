import pandas as pd
##### load config and sample sheets #####
configfile: "config/config.yaml"

sample_table = pd.read_csv(config["sample_table"])
samples = [s for p in pd.concat([sample_table["normal"], sample_table["tumor"]]).str.split(";") for s in p]
patients = [p for p in sample_table["patient"]]

rule s02_01_Mutect2:
    input:
        normal_bam = lambda wildcards: expand(
            "results/01_05_ApplyBQSR/{sample}/{sample}_recal.bam", 
            sample = sample_table.loc[sample_table.patient == wildcards.patient,"normal"]
        ),
        tumor_bam = lambda wildcards: expand(
            "results/01_05_ApplyBQSR/{sample}/{sample}_recal.bam", 
            sample = [s for p in sample_table.loc[sample_table.patient == wildcards.patient,"tumor"].str.split(";") for s in p]
        )   
    output:
        vcf = "results/02_01_Mutect2/{patient}/{patient}_mutect2.vcf.gz",
        f1r2 = "results/02_01_Mutect2/{patient}/{patient}_f1r2.tar.gz",
        stats = "results/02_01_Mutect2/{patient}/{patient}_mutect2.vcf.gz.stats"
    log:
        "results/02_01_Mutect2/{patient}/log" 
    threads: 1
    params:
        ref = config["file_ref_dna_fa"],
        germ_rs = config["germline_resource"], 
        exon_target = config["exon_target"],
        # af = 0.0000025, # --af-of-alleles-not-in-resource {params.af}
        normal = lambda wildcards: expand("{sample}", sample = sample_table.loc[sample_table.patient == wildcards.patient,"normal"]),
        normal_bam = lambda wildcards: expand(
            "-I results/01_05_ApplyBQSR/{sample}/{sample}_recal.bam", 
            sample = sample_table.loc[sample_table.patient == wildcards.patient,"normal"]
        ),
        tumor_bam = lambda wildcards: expand(
            "-I results/01_05_ApplyBQSR/{sample}/{sample}_recal.bam", 
            sample = [s for p in sample_table.loc[sample_table.patient == wildcards.patient,"tumor"].str.split(";") for s in p]
        ) 
    shell:
        '''
        gatk Mutect2 -R {params.ref} {params.tumor_bam} {params.normal_bam} -normal {params.normal} -L {params.exon_target} --germline-resource {params.germ_rs} --f1r2-tar-gz {output.f1r2} -O {output.vcf} 2>{log}
        '''

rule s02_02_LearnReadOrientationModel:
    input:
        rules.s02_01_Mutect2.output.f1r2
    output:
        "results/02_02_LearnReadOrientationModel/{patient}/{patient}_rom.tar.gz"
    log:
        "results/02_02_LearnReadOrientationModel/{patient}/log"
    shell:
        '''
        gatk LearnReadOrientationModel -I {input} -O {output} 2>{log}
        '''

rule s02_03_get_pile:
    input:
        "results/01_05_ApplyBQSR/{sample}/{sample}_recal.bam"
    output:
        "results/02_03_get_pile/{sample}/{sample}_pileup.table"
    log:
        "results/02_03_get_pile/{sample}/log"
    params:
        com_biav = config["common_biallelic_vars"],
        exon_target = config["exon_target"]
    shell:
        '''
        gatk GetPileupSummaries -I {input} -V {params.com_biav} -L {params.exon_target} -O {output} 2>{log}
        '''

rule s02_05_cal_contamination:
    input:
        normal = lambda wildcards: expand(
            "results/02_03_get_pile/{sample}/{sample}_pileup.table", 
            sample = sample_table.loc[sample_table.patient == wildcards.patient,"normal"]
        ),
        tumor = "results/02_03_get_pile/{sample}/{sample}_pileup.table",
    output:
        cont = "results/02_05_cal_contamination/{patient}/{sample}_contamination.table",
        seg = "results/02_05_cal_contamination/{patient}/{sample}_segments.table",
    log:
        "results/02_05_cal_contamination/{patient}/{sample}_log"
    shell:
        '''
        gatk CalculateContamination -I {input.tumor} -matched {input.normal} --tumor-segmentation {output.seg} -O {output.cont} 2>{log}
        '''

rule s02_06_filter:
    input:
        vcf = rules.s02_01_Mutect2.output.vcf,
        f1r2 = rules.s02_01_Mutect2.output.f1r2,
        stats = rules.s02_01_Mutect2.output.stats,
        rom = rules.s02_02_LearnReadOrientationModel.output,
        cont = lambda wildcards: expand(
            "results/02_05_cal_contamination/{patient}/{sample}_contamination.table", 
            patient = wildcards.patient,
            sample = [s for p in sample_table.loc[sample_table.patient == wildcards.patient,"tumor"].str.split(";") for s in p]
        ),
        seg = lambda wildcards: expand(
            "results/02_05_cal_contamination/{patient}/{sample}_segments.table",
            patient = wildcards.patient,
            sample = [s for p in sample_table.loc[sample_table.patient == wildcards.patient,"tumor"].str.split(";") for s in p]
        ) 
    output:
        vcf = "results/02_06_filter/{patient}/{patient}_mutect2_filtered.vcf.gz",
        stats = "results/02_06_filter/{patient}/{patient}_filter.stat"
    log:
        "results/02_06_filter/{patient}/log"
    params:
        ref = config["file_ref_dna_fa"],
        cont = lambda wildcards: expand(
            "--contamination-table results/02_05_cal_contamination/{patient}/{sample}_contamination.table",
            patient = wildcards.patient, 
            sample = [s for p in sample_table.loc[sample_table.patient == wildcards.patient,"tumor"].str.split(";") for s in p]
        ),
        seg = lambda wildcards: expand(
            "--tumor-segmentation results/02_05_cal_contamination/{patient}/{sample}_segments.table",
            patient = wildcards.patient,
            sample = [s for p in sample_table.loc[sample_table.patient == wildcards.patient,"tumor"].str.split(";") for s in p]
        ) 
    shell:
        '''
        gatk FilterMutectCalls -V {input.vcf} -R {params.ref} -O {output.vcf} {params.cont} {params.seg} --ob-priors {input.rom} --stats {input.stats} --filtering-stats {output.stats} 2>{log}
        '''

# rule s02_07_filter_beta:
#     input:
#         vcf = rules.s02_01_Mutect2.output.vcf,
#         f1r2 = rules.s02_01_Mutect2.output.f1r2,
#         stats = rules.s02_01_Mutect2.output.stats,
#         rom = rules.s02_02_LearnReadOrientationModel.output.rom,
#         cont_tb = rules.s02_05_cal_contamination.output.cont_tb,
#         seg_tb = rules.s02_05_cal_contamination.output.seg_tb
#     output:
#         vcf = "results/02_07_filter_beta/{patient}/{tumor}_mutect2_filtered_beta.vcf.gz",
#         stats = "results/02_07_filter_beta/{patient}/{tumor}_filter_beta.stat"
#     log:
#         "results/02_07_filter_beta/{tumor}_VS_{norm}/log"
#     params:
#         ref = config["file_ref_dna_fa"]
#     shell:
#         '''
#         gatk FilterMutectCalls -V {input.vcf} -R {params.ref} -O {output.vcf} --contamination-table {input.cont_tb} --tumor-segmentation {input.seg_tb} --ob-priors {input.rom} --stats {input.stats} --filtering-stats {output.stats} -f-score-beta 1.5 2>{log}
#         '''

rule s02_08_funcotator_vcf:
    input:
        vcf = rules.s02_06_filter.output.vcf
    output:
        vcf = "results/02_08_funcotator_vcf/{patient}/{patient}_mutect2_filtered_funcotator.vcf"
    log:
        "results/02_08_funcotator_vcf/{patient}/log"
    params:
        ref = config["file_ref_dna_fa"],
        exon_target = config["exon_target"],
        db_dir = config["funcotator_db_s"],
        ref_v = config["ref_version"]
    shell:
        '''
        gatk Funcotator -R {params.ref} --ref-version {params.ref_v} -V {input.vcf} -O {output.vcf} -L {params.exon_target} --data-sources-path {params.db_dir} --output-file-format VCF --verbosity DEBUG 2>{log}
        '''

rule s02_09_funcotator_maf:
    input:
        vcf = rules.s02_08_funcotator_vcf.output.vcf
    output:
        maf = "results/02_09_funcotator_maf/{patient}/{patient}_mutect2_filtered_funcotator.maf"
    log:
        "results/02_09_funcotator_maf/{patient}/log"
    params:
        ref = config["file_ref_dna_fa"],
        exon_target = config["exon_target"],
        db_dir = config["funcotator_db_s"],
        ref_v = config["ref_version"]
    shell:
        '''
        gatk Funcotator -R {params.ref} --ref-version {params.ref_v} -V {input.vcf} -O {output.maf} -L {params.exon_target} --data-sources-path {params.db_dir} --output-file-format MAF --verbosity DEBUG 3>{log}
        '''
