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

# genome mapping not cDNA mapping, use dna.primary_assembly.fa not cdna.fa
rule s04_genome_mapping_01_genome_mapping:
    input:
        "results/03_full_length/{sample}/{sample}_full_length.fa",
    output:
        bam = "results/04_genome_mapping/{sample}/{sample}_sorted.bam",
        bai = "results/04_genome_mapping/{sample}/{sample}_sorted.bam.bai",
        done = "results/04_genome_mapping/{sample}/DONE_01_genome_mapping"
    log: "results/04_genome_mapping/{sample}/log_01_genome_mapping"
    threads: 20
    params:
        ref_dna_index = config["file_ref_dna_index"],
        ref_dna_minimap2_bed = config["file_ref_dna_minimap2_bed"],
        fail = "results/04_genome_mapping/{sample}/FAIL"
    shell:
        '''
        minimap2 -t {threads} -ax splice -uf -k14 --secondary=no --junc-bed {params.ref_dna_minimap2_bed} {params.ref_dna_index} {input} | samtools view -q 30 -F 2304 -Sb -|samtools sort - -o {output.bam} 2>{log} && 
        samtools index {output.bam} 2>>{log} &&
        touch {output.done} || touch {params.fail}
        '''

rule s05_quant_known_ref_01_stringtie:
    input:
        "results/04_genome_mapping/{sample}/{sample}_sorted.bam"
    output:
        gene = "results/05_quant_known_ref/{sample}/gene_abund.tab",
        gtf = "results/05_quant_known_ref/{sample}/quant_known_ref_stringtie.gtf",
        done = "results/05_quant_known_ref/{sample}/DONE"
    log: "results/05_quant_known_ref/{sample}/log"
    threads: 4
    params:
        ctab = "results/05_quant_known_ref/{sample}/ctab",
        fail = "results/05_quant_known_ref/{sample}/FAIL",
        gtf = config["file_ref_gtf"],
    shell:
        '''
        stringtie {input} -e -A {output.gene} -b {params.ctab} -p {threads} -G {params.gtf} -o {output.gtf} 2>{log} && 
        touch {output.done} || touch {params.fail}
        '''

rule s05_quant_known_ref_02_stringtie_merge:
    input:
        gtf = expand("results/05_quant_known_ref/{sample}/quant_known_ref_stringtie.gtf", sample = samples),
    output:
        counts = "results/05_quant_known_ref_merge/tran_counts.txt",
    threads: 4
    params:
        sample = samples,
    script:
        "../scripts/s08_quant_corrected_new_ref_02_stringtie_merge.R"

rule s06_build_new_ref_01_merge_part_bam: # avoid excessive files when sambamba merge all bam files
    input:
        bam = lambda wildcards: expand(
            "results/04_genome_mapping/{sample}/{sample}_sorted.bam",
            sample=[s for s in samples if s.startswith(wildcards.part)]
        ),
    output:
        bam = "results/06_build_new_ref/01_part_merged/{part}_merged.bam",
        done = "results/06_build_new_ref/01_part_merged/DONE_01_merge_{part}_bam",
    log: "results/06_build_new_ref/log_01_merge_{part}_bam"
    threads: 20
    shell:
        '''
        sambamba merge -t {threads} -p {output.bam} {input.bam} > {log} 2>&1 &&
        touch {output.done}
        '''

rule s06_build_new_ref_01_merge_all_bam:
    input:
        bam = expand(
            "results/06_build_new_ref/01_part_merged/{part}_merged.bam",
            part = sorted(set(s.split("_")[0] for s in samples))
        ),
    output:
        bam = "results/06_build_new_ref/01_merged.bam",
        done = "results/06_build_new_ref/DONE_01_merge_bam",
    log: "results/06_build_new_ref/log_01_merge_bam"
    threads: 20
    shell:
        '''
        sambamba merge -t {threads} -p {output.bam} {input.bam} > {log} 2>&1 && 
        touch {output.done}
        '''

rule s06_build_new_ref_02_stringtie_assembling:
    input:
        bam = "results/06_build_new_ref/01_merged.bam",
    output:
        gtf = "results/06_build_new_ref/02_stringtie_assembling.gtf",
        done = "results/06_build_new_ref/DONE_02_stringtie_assembling",
    log: "results/06_build_new_ref/log_02_stringtie_assembling"
    threads: 10
    conda:
        "../envs/stringtie.yaml"
    params: 
        gtf = config["file_ref_gtf"],
        name = "CRC"
    shell:
        '''
        stringtie --rf --conservative -L -v -p {threads} -G {params.gtf} -l {params.name} -o {output.gtf} {input.bam} > {log} 2>&1 &&
        touch {output.done} 
        '''

rule s06_build_new_ref_03_sqanti3:
# "cd corresponding working directory" is necessary because the directory of result files is the working directory!!!
    input:
        "results/06_build_new_ref/02_stringtie_assembling.gtf",
    output:
        tmp_gtf = temp("results/06_build_new_ref/temp1.gtf"),
        tmp_cts = temp("results/06_build_new_ref/fl_count.txt"),
        tmp_tx = temp("results/06_build_new_ref/filter_tx.txt"),
        gtf = "results/06_build_new_ref/02_stringtie_assembling_cov3.gtf",
        gff3 = "results/06_build_new_ref/03_sqanti3qc.gff3",
        corrected_fasta = "results/06_build_new_ref/03_sqanti3qc_corrected.fasta",
        corrected_gtf = "results/06_build_new_ref/03_sqanti3qc_corrected.gtf",
        corrected_faa = "results/06_build_new_ref/03_sqanti3qc_corrected.faa",
        classification = "results/06_build_new_ref/03_sqanti3qc_classification.txt",
        filtered_gtf = "results/06_build_new_ref/03_sqanti3filter.filtered.gtf",
        done = "results/06_build_new_ref/DONE_03_sqanti3",
    log:
        "results/06_build_new_ref/log_03_sqanti3"
    threads: 20
    params:
        cdna_cupcake = config["cdna_cupcake"],
        sqanti3_qc = config["dir_sqanti3"] + "/sqanti3_qc.py",
        sqanti3_dir = "results/06_build_new_ref",
        sqanti3_qc_prefix = "03_sqanti3qc",
        cage_peak = config["dir_sqanti3"] + "/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed",
        polya_motif= config["dir_sqanti3"] + "/data/polyA_motifs/mouse_and_human.polyA_motif.txt",
        tappas_gff = "resources/tappAS/Homo_sapiens_GRCh38_Ensembl_86.gff3",
        ref_gtf = config["file_ref_gtf"],
        ref_fa = config["file_ref_dna_fa"],
        sqanti3_filter = config["dir_sqanti3"] + "/sqanti3_filter.py",
        sqanti3_filter_prefix = "03_sqanti3filter",
    shell:
        '''
        ########### filtering isoforms without strand information and coverage <= 3
        cat {input} | awk '$7!="." {{print $0}}' > {output.tmp_gtf} &&
        echo "pbid\tcount_fl" > {output.tmp_cts} &&
        cat {output.tmp_gtf} |grep -w transcript |grep ref_gene_name |sed 's/\"\|\;//g'|cut -d ' ' -f 4,12 |awk '{{print $1,int($2)}}' | awk '{{if($2>3) print $0}}'| sed 's/\s/\\t/g' >> {output.tmp_cts} &&
        cat {output.tmp_gtf} |grep -w transcript |grep -v ref_gene_name |sed 's/\"\|\;//g'|cut -d ' ' -f 4,6 |awk '{{print $1,int($2)}}' | awk '{{if($2>3) print $0}}'| sed 's/\s/\\t/g' >> {output.tmp_cts} &&
        tail -n +2 {output.tmp_cts} | cut -f 1 > {output.tmp_tx} &&
        grep -F -f {output.tmp_tx} {output.tmp_gtf} > {output.gtf} &&
        ########### sqanti3 qc
        {params.cdna_cupcake}
        python {params.sqanti3_qc} -t {threads} -d {params.sqanti3_dir} -o {params.sqanti3_qc_prefix} --report pdf --CAGE_peak {params.cage_peak} --polyA_motif_list {params.polya_motif} --isoAnnotLite --gff3 {params.tappas_gff} {output.gtf} {params.ref_gtf} {params.ref_fa} > {log} 2>&1 &&
        ########### sqanti3 filter
        python {params.sqanti3_filter} ML --isoAnnotGFF3 {output.gff3} --isoforms {output.corrected_fasta} --gtf {output.corrected_gtf} --faa {output.corrected_faa} -d {params.sqanti3_dir} -o {params.sqanti3_filter_prefix} {output.classification} >> {log} 2>&1 &&
        touch {output.done}
        '''

rule s07_quant_new_ref_01_stringtie:
    input:
        bam = "results/04_genome_mapping/{sample}/{sample}_sorted.bam",
        gtf = "results/06_build_new_ref/03_sqanti3filter.filtered.gtf",
    output:
        gene = "results/07_quant_new_ref/{sample}/gene_abund.tab",
        gtf = "results/07_quant_new_ref/{sample}/quant_new_ref_stringtie.gtf",
        done = "results/07_quant_new_ref/{sample}/DONE"
    log: "results/07_quant_new_ref/{sample}/log"
    threads: 4
    params:
        ctab = "results/07_quant_new_ref/{sample}/ctab",
        fail = "results/07_quant_new_ref/{sample}/FAIL"
    shell:
        '''
        stringtie {input.bam} -e -A {output.gene} -b {params.ctab} -p {threads} -G {input.gtf} -o {output.gtf} 2>{log} && 
        touch {output.done} || touch {params.fail}
        '''

rule s07_quant_new_ref_02_stringtie_merge:
    input:
        gtf = expand("results/07_quant_new_ref/{sample}/quant_new_ref_stringtie.gtf", sample = samples)
    output:
        counts = "results/07_quant_new_ref_merge/tran_counts.txt",
    threads: 4
    params:
        sample = samples,
    script:
        "../scripts/s07_quant_new_ref_02_stringtie_merge.R"

rule s08_quant_corrected_new_ref_01_stringtie:
    input:
        bam = "results/04_genome_mapping/{sample}/{sample}_sorted.bam",
        gtf = "results/08_corrected_new_ref/corrected_new_ref.sorted.gtf", # scripts/08_corrected_new_ref.ipynb
    output:
        gene = "results/08_quant_corrected_new_ref/{sample}/gene_abund.tab",
        gtf = "results/08_quant_corrected_new_ref/{sample}/quant_corrected_new_ref_stringtie.gtf",
        done = "results/08_quant_corrected_new_ref/{sample}/DONE"
    log: "results/08_quant_corrected_new_ref/{sample}/log"
    threads: 4
    params:
        ctab = "results/08_quant_corrected_new_ref/{sample}/ctab",
        fail = "results/08_quant_corrected_new_ref/{sample}/FAIL"
    shell:
        '''
        stringtie {input.bam} -e -A {output.gene} -b {params.ctab} -p {threads} -G {input.gtf} -o {output.gtf} 2>{log} && 
        touch {output.done} || touch {params.fail}
        '''

rule s08_quant_corrected_new_ref_02_stringtie_merge:
    input:
        gtf = expand("results/08_quant_corrected_new_ref/{sample}/quant_corrected_new_ref_stringtie.gtf", sample = samples)
    output:
        counts = "results/08_quant_corrected_new_ref_merge/tran_counts.txt",
    threads: 4
    params:
        sample = samples,
    script:
        "../scripts/s08_quant_corrected_new_ref_02_stringtie_merge.R"
