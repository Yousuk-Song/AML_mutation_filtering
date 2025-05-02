#!/usr/bin/bash

# 참조 유전체 및 툴 경로
ref=/data/processed_data/scRSEQ_AML/MUTECT_MUTATION/ref/GRCh38.p13.genome.fa
gatk=/data/program/gatk/gatk-4.3.0.0/gatk

# 입력 BAM 파일
bam=$1

# Panel of Normals (PON) 및 Germline Resource 경로
pon=/data/processed_data/scRSEQ_AML/MUTECT_MUTATION/ref/1000g_pon.hg38.vcf.gz
germline_resource=/data/processed_data/scRSEQ_AML/MUTECT_MUTATION/ref/af-only-gnomad.hg38.vcf.gz

# Mutect2 인터벌 (bed)
pancancer=/data/processed_data/scRSEQ_AML/marker_macrogen/Adjusted_bam_dir/Pancancer/pan_cancer_coordinates_sorted.bed

# @RG 태그에서 SM 추출
tumor_sample=$(samtools view -H ${bam} | grep '@RG' | sed -n 's/.*SM:\([^\t]*\).*/\1/p')
echo "Tumor Sample Name: ${tumor_sample}"

if [[ -z "${tumor_sample}" ]]; then
    echo "Error: Unable to extract tumor sample name from BAM file header"
    exit 1
fi

# Mutect2 실행
$gatk Mutect2 \
   -R ${ref} \
   -I ${bam} \
   -tumor ${tumor_sample} \
   --germline-resource ${germline_resource} \
   --panel-of-normals ${pon} \
   --intervals ${pancancer} \
   -O ${bam}.mutect2_pancancer.vcf
