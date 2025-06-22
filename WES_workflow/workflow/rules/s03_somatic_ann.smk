import pandas as pd
##### load config and sample sheets #####
configfile: "config/config.yaml"

sample_table = pd.read_csv(config["sample_table"])
samples = [s for p in pd.concat([sample_table["normal"], sample_table["tumor"]]).str.split(";") for s in p]
patients = [p for p in sample_table["patient"]]


rule s03_01_snpEff:
    input:
        vcf = "results/02_06_filter/{patient}/{patient}_mutect2_filtered.vcf.gz"
    output:
        vcf = "results/03_01_snpEff/{patient}/{patient}_mutect2_filtered_snpEff.vcf"
    log:
        "results/03_01_snpEff/{patient}/log"
    params:
        snpEff = config["SOFTWARE"]["snpEff"],
        snpEff_config = config["SOFTWARE"]["snpEff_config"]
    shell:
        '''
        java -Xmx4g -jar {params.snpEff} -c {params.snpEff_config} hg38 {input.vcf} > {output.vcf} 2>{log}
        '''

rule s03_02_snpSift:
    input:
        vcf = rules.s03_01_snpEff.output
    output:
        vcf = "results/03_02_snpSift/{patient}/{patient}_mutect2_filtered_snpSift.vcf"
    log:
        "results/03_02_snpSift/{patient}/log"
    params:
        snpSift = config["SOFTWARE"]["snpSift"],
        snpSift_db = config["SOFTWARE"]["snpSift_db"]
    shell:
        '''
        java -Xmx4g -jar {params.snpSift} dbnsfp -v -db {params.snpSift_db} {input.vcf} > {output.vcf} 2>{log}
        '''

rule s03_03_left_norm:
    input:
        vcf = "results/02_06_filter/{patient}/{patient}_mutect2_filtered.vcf.gz"
    output:
        vcf = "results/03_03_left_norm/{patient}/{patient}_mutect2_filtered_av.vcf"
    log:
        "results/03_03_left_norm/{patient}/log"
    params:
        ref_fa = config["file_ref_dna_fa"]
    shell:
        '''
        bcftools norm -m -any {input} | bcftools norm -f {params.ref_fa} -o {output.vcf} -
        '''

rule s03_04_annovar:
    input:
        vcf = rules.s03_03_left_norm.output.vcf
    output:
        vcf = "results/03_04_annovar/{patient}/{patient}.hg38_multianno.vcf"
    log:
        "results/03_04_annovar/{patient}/log"
    params:
        exe_annovar = config["SOFTWARE"]["annovar"],
        ann_database = config["SOFTWARE"]["annovar_database"],
        out = "results/03_04_annovar/{patient}/{patient}",
        gene_ann = config["SOFTWARE"]["annovar_gene_ann"] 
    shell:
        '''
        {params.exe_annovar} {input.vcf} {params.ann_database} -buildver hg38 -out {params.out} -remove -protocol refGene,cytoBand,exac03,avsnp150,dbnsfp41a -operation gx,r,f,f,f -nastring . -polish -vcfinput 2>{log}
        '''

rule s03_04_annovar_vcf2maf:
    input:
        "results/03_04_annovar/{patient}/{patient}.hg38_multianno.vcf"
    output:
        "results/03_04_annovar/{patient}/{sample}/{sample}.maf"
    log:
        "results/03_04_annovar/{patient}/{sample}/{sample}.log"
    params:
        tmp_dir = "results/03_04_annovar/{patient}/{sample}",
        tumor_id = "{sample}",
        norm_id = lambda wildcards: expand("{sample}", sample = sample_table.loc[sample_table.patient == wildcards.patient,"normal"]),
    shell:
        '''
        vcf2maf.pl --input-vcf {input} --output-maf {output} --tmp-dir {params.tmp_dir} --tumor-id {params.tumor_id} --normal-id {params.norm_id} --vep-data /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/workflow/exon_snv/database/vep --vep-path /gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/luping/software/conda/envs/gatk/bin/ --verbose --ref-fasta ~/workflow/exon_snv/database/GRCh38/GRCh38_chr.dna.fa --ncbi-build GRCh38
        ''' 


