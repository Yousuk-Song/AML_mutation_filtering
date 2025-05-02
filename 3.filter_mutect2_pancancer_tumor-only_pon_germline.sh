#!/usr/bin/bash

ref=/data/processed_data/scRSEQ_AML/MUTECT_MUTATION/ref/GRCh38.p13.genome.fa
gatk=/data/program/gatk/gatk-4.3.0.0/gatk


INPUT_VCF=$1
CONTABLE=$2

OUTPUT_VCF="$(basename "$INPUT_VCF" .vcf).filtered.vcf"

$gatk FilterMutectCalls \
   -V ${INPUT_VCF} \
   -R ${ref} \
   --contamination-table ${CONTABLE} \
   --max-events-in-region 5 \
   -O ${OUTPUT_VCF} 

echo "Filtered VCF file generated: ${OUTPUT_VCF}"
