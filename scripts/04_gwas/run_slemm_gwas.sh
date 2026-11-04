#!/usr/bin/env bash
set -euo pipefail
: "${GENO_DIR:?Set the authorized genotype prefix}"
: "${PHENO_FILE:?Set the phenotype file}"
: "${OUTPUT_DIR:?Set a NEW output directory}"
THREADS=${THREADS:-16}
[[ ! -e "$OUTPUT_DIR" ]] || { echo 'Output directory must be new' >&2; exit 2; }
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
: "${SNP_INFO_FILE:?Set the historical SLEMM SNP information file}"
SLEMM_BIN=${SLEMM_BIN:-slemm}
SLEMM_GAMMA_BIN=${SLEMM_GAMMA_BIN:-slemm_gamma}
SLEMM_GWA_BIN=${SLEMM_GWA_BIN:-slemm_gwa}
export OMP_NUM_THREADS="$THREADS"
for i in {1..10}; do
  PREFIX="$OUTPUT_DIR/spg.even.100k.pheno$i"
  "$SLEMM_BIN" --lmm --phenotype_file "$PHENO_FILE" --trait "pheno$i"     --bfile "$GENO_DIR" --snp_info_file "$SNP_INFO_FILE" --out "$PREFIX"
  OMP_NUM_THREADS=1 "$SLEMM_GAMMA_BIN" --pfile "$GENO_DIR" --slemm "$PREFIX" --out "$PREFIX.gamma.txt"
  for chr in {1..18}; do
    "$SLEMM_GWA_BIN" --pfile "$GENO_DIR" --slemm "$PREFIX" --chr "$chr" --out "$PREFIX.chr$chr.gwas.txt"
  done
  cp "$PREFIX.chr1.gwas.txt" "$PREFIX.allchr.gwas.txt"
  for chr in {2..18}; do tail -n +2 "$PREFIX.chr$chr.gwas.txt" >> "$PREFIX.allchr.gwas.txt"; done
done
# Keep per-chromosome outputs for auditing. No cleanup of originals.
