#!/usr/bin/env bash
set -euo pipefail
: "${GENO_DIR:?Set the authorized genotype prefix}"
: "${PHENO_FILE:?Set the phenotype file}"
: "${OUTPUT_DIR:?Set a NEW output directory}"
THREADS=${THREADS:-16}
[[ ! -e "$OUTPUT_DIR" ]] || { echo 'Output directory must be new' >&2; exit 2; }
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
REGENIE_BIN=${REGENIE_BIN:-regenie}
BSIZE=${BSIZE:-1000}
for i in {1..10}; do
  "$REGENIE_BIN" --step 1 --bed "$GENO_DIR" --phenoFile "$PHENO_FILE" --phenoCol "pheno$i"     --qt --bsize "$BSIZE" --threads "$THREADS" --out "$OUTPUT_DIR/regenie_step1_pheno$i"
  "$REGENIE_BIN" --step 2 --bed "$GENO_DIR" --phenoFile "$PHENO_FILE" --phenoCol "pheno$i"     --qt --bsize "$BSIZE" --threads "$THREADS" --pred "$OUTPUT_DIR/regenie_step1_pheno${i}_pred.list"     --out "$OUTPUT_DIR/regenie_pheno$i"
done
