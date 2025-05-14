#!/usr/bin/env python3

import pysam
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import sys
import os

# === 파일 경로 지정 ===
tumor_vcf_path = sys.argv[1]
normal_vcf_path = sys.argv[2]

# === 샘플명 추출 (ex. pt01_bm01)
sample_prefix = os.path.basename(tumor_vcf_path).split("_tumor")[0]
hist_output = f"{sample_prefix}_tumor_normal_AF_difference_histogram.png"
venn_output = f"{sample_prefix}_tumor_normal_AF_venn.png"

# === AF 추출 함수 ===
def extract_af_table(vcf_path):
    vcf = pysam.VariantFile(vcf_path)
    records = []
    for rec in vcf.fetch():
        chrom = rec.chrom
        pos = rec.pos
        ref = rec.ref
        alt = rec.alts[0] if rec.alts else None
        af = rec.samples[0].get('AF', [None])[0]
        if af is not None:
            records.append((chrom, pos, ref, alt, af))
    return pd.DataFrame(records, columns=["CHROM", "POS", "REF", "ALT", "AF"])

# === 데이터 로딩 ===
tumor_df = extract_af_table(tumor_vcf_path).rename(columns={"AF": "AF_tumor"})
normal_df = extract_af_table(normal_vcf_path).rename(columns={"AF": "AF_normal"})

# === 변이 매칭 및 차이 계산 ===
merged = pd.merge(tumor_df, normal_df, on=["CHROM", "POS", "REF", "ALT"])
merged["AF_diff"] = merged["AF_tumor"] - merged["AF_normal"]

# === 히스토그램 출력 ===
plt.figure(figsize=(10, 6))
plt.hist(merged["AF_diff"], bins=30, color="purple", edgecolor="black", alpha=0.8)
plt.axvline(x=0, linestyle="--", color="gray")
plt.title(f"{sample_prefix}: AF Difference (Tumor - Normal)")
plt.xlabel("AF Difference")
plt.ylabel("Number of Shared Mutations")
plt.tight_layout()
plt.savefig(hist_output, dpi=150)
plt.close()

# === Venn 다이어그램 출력 ===
tumor_set = set(zip(tumor_df["CHROM"], tumor_df["POS"], tumor_df["REF"], tumor_df["ALT"]))
normal_set = set(zip(normal_df["CHROM"], normal_df["POS"], normal_df["REF"], normal_df["ALT"]))

plt.figure(figsize=(6, 6))
venn2([tumor_set, normal_set], set_labels=("Tumor", "Normal"))
plt.title(f"{sample_prefix}: Variant Overlap (CHROM,POS,REF,ALT)")
plt.tight_layout()
plt.savefig(venn_output, dpi=150)
plt.close()
