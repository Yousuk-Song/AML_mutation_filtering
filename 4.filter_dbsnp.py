#!/usr/bin/python

import pysam
import os
import sys

if len(sys.argv) < 2:
    print("❗ 사용법: python filter_dbsnp.py <input.vcf>")
    sys.exit(1)

# 입력 파일
input_vcf_path = sys.argv[1]
dbsnp_vcf_path = "/data/processed_data/scRSEQ_AML/MUTECT_MUTATION/ref/GCF_000001405.40.af_ge_0.5.cleaned.vcf.gz"
output_vcf_path = input_vcf_path.replace(".vcf", ".dbsnp_filtered.vcf")

# dbSNP 인덱스 존재 확인
if not os.path.exists(dbsnp_vcf_path + ".tbi"):
    raise FileNotFoundError("❗ dbSNP VCF는 반드시 bgzip + index 되어 있어야 합니다 (.tbi 필요).")

# 파일 열기

input_vcf = pysam.VariantFile(input_vcf_path)
dbsnp_vcf = pysam.TabixFile(dbsnp_vcf_path)
output_vcf = pysam.VariantFile(output_vcf_path, "w", header=input_vcf.header)

# 필터링 수행
for record in input_vcf.fetch():
    chrom = record.chrom
    pos = record.pos
    ref = record.ref
    alts = record.alts  # tuple, 예: ("C", "G")

    keep = True
    try:
        for db_line in dbsnp_vcf.fetch(chrom, pos - 1, pos):
            fields = db_line.strip().split("\t")
            db_ref = fields[3]
            db_alts = fields[4].split(",")

            if ref == db_ref:
                if any(alt in db_alts for alt in alts):
                    keep = False
                    break
    except Exception:
        pass  # 해당 구간에 dbSNP 없을 경우 무시

    if keep:
        output_vcf.write(record)

# 마무리
input_vcf.close()
output_vcf.close()

# bgzip 압축
subprocess.run(["bgzip", "-f", output_vcf_path], check=True)

# tbi 인덱스 생성
subprocess.run(["tabix", "-p", "vcf", output_vcf_gz], check=True)

print(f"✅ 필터링 및 압축 완료: {output_vcf_gz}")
print(f"✅ 인덱스 생성 완료: {output_vcf_gz}.tbi")
