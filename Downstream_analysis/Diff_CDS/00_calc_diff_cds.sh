# conda activate pygenometracks
# step 1
awk '{if ($3 == "CDS") {print $0}}' ~/project/ont_iso_crc/2_RNA_CRC/0_smk_CRC_HTCRC/results/08_corrected_new_ref/CRC_HTCRC_corrected_new_ref.sorted.gtf | sort -k1,1 -k4,4n > 1_cds.gtf
sort -k1,1 -k4,4n ~/project/ont_iso_crc/2_RNA_CRC/0_smk_CRC_HTCRC/resources/hg38_Pfam_in_GENCODE.gtf | sed 's/^chr//g' > 1_pfam.gtf

# step 2
python 2_cds_diff.py -p 01_DTU_tx_each_patient_AllInfo_CRC01.txt -g 1_cds.gtf -o 2_cds_diff.bed -j8
sort -k1,1 -k2,2n 2_cds_diff.bed > 2_cds_diff.sorted.bed

# step 3
bedtools intersect -a 2_cds_diff.sorted.bed -b 1_pfam.gtf -wb > 3_cds_diff_in_pfam.txt

# step 4
python 4_summarize_cds_diff.py  -i 2_cds_diff.bed -o 4_cds_diff_summary.txt
python 4_summarize_pfam_diff.py -i 3_cds_diff_in_pfam.txt -o 4_cds_diff_in_pfam_summary.txt
