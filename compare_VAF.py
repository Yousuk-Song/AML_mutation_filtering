#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
vcf_normal <- args[1]
vcf_tumor <- args[2]

library(VariantAnnotation)
library(ggplot2)

# === 파일명, 라벨 설정 ===
extract_base_label <- function(filename) {
  sub("\\.withSM.*", "", filename)
}

sample_name_normal <- extract_base_label(basename(vcf_normal))
sample_name_tumor  <- extract_base_label(basename(vcf_tumor))

# 파일명에 pathogenic 들어 있으면 flag 붙이기
is_pathogenic <- grepl("pathogenic", vcf_normal, ignore.case = TRUE) |
                 grepl("pathogenic", vcf_tumor,  ignore.case = TRUE)

file_suffix <- if (is_pathogenic) {
  ".pathogenic.vaf_histogram_ggplot.png"
} else {
  ".vaf_histogram_ggplot.png"
}

# 플롯 제목
plot_title <- paste0(sample_name_normal, " & ", sample_name_tumor,
                     if (is_pathogenic) " (pathogenic)" else "")

output_file <- paste0(sample_name_normal, "_vs_", sample_name_tumor, file_suffix)

# === VAF 추출 함수 ===
extract_vaf <- function(vcf_file) {
  vcf <- readVcf(vcf_file, genome = "hg38")
  vcf_valid <- vcf[fixed(vcf)$FILTER == "PASS" | grepl("pathogenic", fixed(vcf)$FILTER, ignore.case = TRUE)]
  vafs <- as.numeric(sapply(geno(vcf_valid)$AF, function(x) x[1]))
  vafs[!is.na(vafs) & vafs >= 0 & vafs <= 1]
}

# === 데이터 불러오기 ===
vaf_normal <- extract_vaf(vcf_normal)
vaf_tumor  <- extract_vaf(vcf_tumor)

df <- data.frame(
  VAF = c(vaf_normal, vaf_tumor),
  Group = factor(c(rep(sample_name_normal, length(vaf_normal)),
                   rep(sample_name_tumor, length(vaf_tumor))))
)

# === 중앙값 계산 ===
vaf_medians <- aggregate(VAF ~ Group, data = df, FUN = median)

# === 그룹별 히스토그램 count 최대값 ===
ymax_count_normal <- max(hist(vaf_normal, breaks = 30, plot = FALSE)$counts)
ymax_count_tumor  <- max(hist(vaf_tumor,  breaks = 30, plot = FALSE)$counts)

# === 그룹별 density 계산 및 개별 스케일링 ===
dens_normal <- density(vaf_normal)
dens_tumor  <- density(vaf_tumor)

scale_normal <- ymax_count_normal / max(dens_normal$y)
scale_tumor  <- ymax_count_tumor / max(dens_tumor$y)

df_dens <- rbind(
  data.frame(VAF = dens_normal$x, y = dens_normal$y * scale_normal, Group = sample_name_normal),
  data.frame(VAF = dens_tumor$x,  y = dens_tumor$y  * scale_tumor,  Group = sample_name_tumor)
)

# === 중앙값 텍스트 y 좌표 설정 ===
ymax_total <- max(c(ymax_count_normal, ymax_count_tumor))
vaf_medians$y <- c(ymax_total * 0.9, ymax_total * 0.75)
vaf_medians$Group <- factor(vaf_medians$Group)

# === 플롯 ===
p <- ggplot(df, aes(x = VAF, fill = Group)) +
  geom_histogram(aes(y = after_stat(count), color = Group),
                 bins = 30, alpha = 0.35, position = "identity") +
  geom_line(data = df_dens, aes(x = VAF, y = y, color = Group), linewidth = 1.2) +
  geom_vline(data = vaf_medians, aes(xintercept = VAF, color = Group), linetype = "dashed") +
  geom_text(data = vaf_medians,
            aes(x = VAF, y = y, label = paste0("Median: ", round(VAF, 2)), color = Group),
            inherit.aes = FALSE, hjust = -0.1, size = 5) +
  scale_fill_manual(values = c("dodgerblue", "red")) +
  scale_color_manual(values = c("dodgerblue", "red")) +
  scale_y_continuous(
    name = "Mutation Count",
    sec.axis = sec_axis(~ ., name = "Scaled Density (per group)")
  ) +
  labs(
    title = plot_title,
    x = "VAF"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    axis.title.y.right = element_text(color = "black")
  )

# === 저장 ===
ggsave(output_file, plot = p, width = 10, height = 8, dpi = 150)
