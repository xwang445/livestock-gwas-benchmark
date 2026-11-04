#!/usr/bin/env bash
set -euo pipefail
: "${GENO_DIR:?Set the authorized genotype prefix}"
: "${PHENO_FILE:?Set the phenotype file}"
: "${OUTPUT_DIR:?Set a NEW output directory}"
THREADS=${THREADS:-16}
[[ ! -e "$OUTPUT_DIR" ]] || { echo 'Output directory must be new' >&2; exit 2; }
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
BOLT_BIN=${BOLT_BIN:-bolt}
for i in {1..10}; do
  "$BOLT_BIN" --bfile "$GENO_DIR" --phenoFile "$PHENO_FILE" --phenoCol "pheno$i"     --lmm --LDscoresUseChip --numThreads "$THREADS" --statsFile "$OUTPUT_DIR/bolt_pheno$i.stats"
done
