#!/usr/bin/env Rscript

# === 입력 처리 ===
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
output_file <- paste0(vcf_file, ".vaf_histogram_ggplot.png")

# === 라이브러리 로드 ===
library(VariantAnnotation)
library(ggplot2)

# === 데이터 처리 ===
vcf <- readVcf(vcf_file, genome = "hg38")
vcf_valid <- vcf[fixed(vcf)$FILTER == "PASS" | grepl("pathogenic", fixed(vcf)$FILTER, ignore.case = TRUE)]
vafs <- as.numeric(sapply(geno(vcf_valid)$AF, function(x) x[1]))
vafs <- vafs[!is.na(vafs) & vafs >= 0 & vafs <= 1]

df <- data.frame(VAF = vafs)
sample_name <- gsub("_Bioproduct.*", "", basename(vcf_file))
vaf_median <- median(df$VAF)

# === 그래프 생성 ===
p <- ggplot(df, aes(x = VAF)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 20, 
                 fill = "dodgerblue", 
                 color = "white", 
                 alpha = 0.8) +
  geom_density(color = "red", linewidth = 1.2, adjust = 1) +
  geom_vline(xintercept = vaf_median, 
             linetype = "dashed", 
             color = "darkgreen", 
             linewidth = 1) +
  annotate("text", 
           x = vaf_median, 
           y = max(density(df$VAF)$y) * 0.95, 
           label = paste0("Median: ", round(vaf_median, 2)),
           color = "darkgreen",
           size = 5,
           hjust = -0.1) +
  labs(
    title = "VAF Histogram with Density",
    subtitle = sample_name,
    x = "VAF",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# === 저장 ===
ggsave(output_file, plot = p, width = 10, height = 8, dpi = 150)
