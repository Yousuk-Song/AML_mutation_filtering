#!/usr/bin/env Rscript

# 명령줄 인자에서 VCF 파일 경로를 받아옴
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]

# 출력 PNG 파일 이름 정의
output_file <- paste0(vcf_file, ".vaf_density.png")

# 필요한 패키지 로드
library(VariantAnnotation)

# VCF 파일 로드 (압축된 .vcf.gz 가능)
vcf <- readVcf(vcf_file, genome = "hg38")

# PASS된 변이만 필터링
vcf_pass <- vcf[fixed(vcf)$FILTER == "PASS"]

# VAF 추출: 형식 필드 FORMAT에서 AF 값을 numeric으로 추출
vafs <- as.numeric(sapply(geno(vcf_pass)$AF, function(x) x[1]))

# VAF 값이 NA가 아닌 것만 사용
vafs <- vafs[!is.na(vafs)]

# PNG 파일로 저장 (type = "cairo" 설정으로 X11 없는 서버에서 사용 가능)
png(filename = output_file, width = 800, height = 600, type = "cairo")

# VAF density plot 생성
plot(density(vafs), main = "VAF Density (PASS variants)", xlab = "VAF", ylab = "Density")

# 장치 종료 (파일 저장 완료)
dev.off()
