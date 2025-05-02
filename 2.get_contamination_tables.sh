#!/usr/bin/bash

# 환경 변수 설정
VCF="/data/processed_data/StepwiseHCC_WGS/Mutec2/resource/small_exac_common_3.hg38.vcf.gz"
INTERVALS="/data/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/Pancancer/pan_cancer_coordinates_sorted.bed"
gatk=/data/program/gatk/gatk-4.3.0.0/gatk

# BAM 파일 루프
bam=$1

echo "▶️ Processing $sample..."

# Step 1: GetPileupSummaries
$gatk GetPileupSummaries \
-I "$bam" \
-V "$VCF" \
3-L "$INTERVALS" \
-O contamination_table_dir/"${bam}.pileups.table" 

# Step 2: CalculateContamination
$gatk CalculateContamination \
-I contamination_table_dir/"${bam}.pileups.table" \
-O contamination_table_dir/"${bam}.contamination.table"
