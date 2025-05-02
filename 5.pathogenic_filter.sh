#!/usr/bin/env bash

# --- Required tools ---
# 필요한 도구: bcftools, bgzip, tabix, annovar (convert2annovar.pl, table_annovar.pl)

# --- Input ---
input_vcf=$1  # 첫 번째 인자로 입력된 VCF.gz 파일 경로 (예: sample.vcf.gz)
annovar_db="/data/program/annovar/humandb"  # ANNOVAR 데이터베이스 디렉토리
ref_ver="hg38"  # 참조 유전체 버전
threads=4  # 멀티스레드 수 설정

# --- Check ---
if [[ ! -f "$input_vcf" ]]; then
    echo "❌ VCF not found: $input_vcf"
    exit 1  # 입력 파일 없으면 종료
fi

# --- Derived filenames ---
base=$(basename "$input_vcf" .vcf.gz)  # 입력 파일명에서 ".vcf.gz" 제거한 베이스 이름 추출
out_dir="./pathogenic_filter_work"  # 결과 저장 디렉토리
mkdir -p "$out_dir"  # 결과 디렉토리 생성 (존재하지 않으면)

filtered_vcf="${out_dir}/${base}.pass_or_germline.vcf.gz"  # PASS 또는 germline 변이만 포함할 VCF 경로
annovar_input="${out_dir}/${base}.avinput"  # ANNOVAR 입력 포맷 (.avinput) 파일 경로
multianno_csv="${out_dir}/${base}.${ref_ver}_multianno.csv"  # ANNOVAR 결과 CSV 파일 경로
patho_ids="${out_dir}/${base}.pathogenic.ids.txt"  # 병리적 변이의 chr 및 pos만 추출한 ID 리스트
patho_vcf="${out_dir}/${base}.pathogenic.vcf.gz"  # 병리적 변이만 포함한 최종 VCF 경로

# [1/5] PASS 또는 germline 필터 통과 변이만 추출
echo "[1/5] Extracting PASS/germline variants..."
bcftools view -f PASS,germline "$input_vcf" -Oz -o "$filtered_vcf" --threads $threads
tabix -p vcf "$filtered_vcf"  # 인덱스 생성

# [2/5] VCF → ANNOVAR 포맷(.avinput)으로 변환
echo "[2/5] Converting to ANNOVAR input..."
convert2annovar.pl -format vcf4old "$filtered_vcf" > "$annovar_input"

# [3/5] ANNOVAR 기능 유전자 및 ClinVar 정보로 주석 부여
echo "[3/5] Running ANNOVAR annotation..."
table_annovar.pl "$annovar_input" "$annovar_db" \
  -buildver "$ref_ver" \
  -out "${out_dir}/${base}" \
  -remove \  # 임시 파일 삭제
  -protocol refGene,clinvar_20220320 \  # 유전자 및 ClinVar 병리성 주석 사용
  -operation g,f \  # g: gene-based, f: filter-based
  -nastring . \  # 결측치는 '.'으로 표기
  -csvout  # CSV 형식으로 출력

# [4/5] CSV에서 "Pathogenic" 또는 "Likely_pathogenic" 라벨이 붙은 변이만 추출
echo "[4/5] Extracting pathogenic variants from annotated CSV..."
awk -F',' 'NR==1 || $0 ~ /Pathogenic|Likely_pathogenic/' "$multianno_csv" | \
    awk -F',' 'NR>1 {print $1"\t"$2}' > "$patho_ids"
# 첫 줄(헤더) 제외 후, chr와 pos 추출 → 이후 VCF 추출용 위치 리스트 생성

# [5/5] 해당 위치(patho_ids)에 해당하는 변이만 최종 VCF에서 추출
echo "[5/5] Extracting pathogenic variants from VCF..."
bcftools view -R "$patho_ids" "$filtered_vcf" -Oz -o "$patho_vcf" --threads $threads
tabix -p vcf "$patho_vcf"  # 최종 VCF 인덱스 생성

echo "✅ Final pathogenic-only VCF: $patho_vcf"
