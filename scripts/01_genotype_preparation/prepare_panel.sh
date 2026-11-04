#!/usr/bin/env bash
set -euo pipefail
# Reconstruct the observed QC + thinning + QC-passing-chip union.
# Source pgens, CHIP_IDS, and historical SNP/sample order are restricted inputs.
: "${SOURCE_PGEN_DIR:?Set the directory containing chr1..18.swim.imputed pgens}"
: "${CHIP_IDS:?Set a one-ID-per-line list of chip SNPs; not the causal-only list}"
: "${OUTPUT_DIR:?Set a NEW output directory}"
PLINK2_BIN=${PLINK2_BIN:-plink2}
THREADS=${THREADS:-20}
[[ ! -e "$OUTPUT_DIR" ]] || { echo 'Output must not already exist' >&2; exit 2; }
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
for chr in {1..18}; do
  "$PLINK2_BIN" --pfile "$SOURCE_PGEN_DIR/chr${chr}.swim.imputed" \
    --extract-if-info 'INFO >= 0.8' --maf 0.01 --hwe 1e-9 \
    --make-pgen --write-snplist --threads "$THREADS" --out "$OUTPUT_DIR/chr${chr}.qc"
  "$PLINK2_BIN" --pfile "$OUTPUT_DIR/chr${chr}.qc" --bp-space 20000 \
    --write-snplist --threads "$THREADS" --out "$OUTPUT_DIR/chr${chr}.even"
done
cat "$OUTPUT_DIR"/chr{1..18}.even.snplist > "$OUTPUT_DIR/all_chr_even.snplist"
# Preserve first-occurrence order like unique(c(even_snps, chip_snps)) in R.
awk '!seen[$0]++' "$OUTPUT_DIR/all_chr_even.snplist" "$CHIP_IDS" > "$OUTPUT_DIR/snp.final.snplist"
for chr in {1..18}; do printf '%s\n' "$OUTPUT_DIR/chr${chr}.qc"; done > "$OUTPUT_DIR/merge_list.txt"
"$PLINK2_BIN" --pmerge-list "$OUTPUT_DIR/merge_list.txt" --extract "$OUTPUT_DIR/snp.final.snplist" \
  --make-pgen --threads "$THREADS" --out "$OUTPUT_DIR/spg.even.100k"
"$PLINK2_BIN" --pfile "$OUTPUT_DIR/spg.even.100k" --make-bed --out "$OUTPUT_DIR/spg.even.100k"
[[ $(wc -l < "$OUTPUT_DIR/spg.even.100k.bim") -eq 135358 ]]
[[ $(wc -l < "$OUTPUT_DIR/spg.even.100k.fam") -eq 31301 ]]
# Intermediate files are retained. Exact original chip-list contents are missing;
# the recovered union rule reproduces observed membership, not byte identity.
