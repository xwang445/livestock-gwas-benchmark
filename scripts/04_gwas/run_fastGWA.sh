#!/usr/bin/env bash
set -euo pipefail
: "${GENO_DIR:?Set the authorized genotype prefix}"
: "${PHENO_DIR:?Set fastGWA_data directory}"
: "${OUTPUT_DIR:?Set a NEW output directory}"
THREADS=${THREADS:-16}
[[ ! -e "$OUTPUT_DIR" ]] || { echo 'Output directory must be new' >&2; exit 2; }
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
GCTA_BIN=${GCTA_BIN:-gcta_1.94.1}
PREFIX="$OUTPUT_DIR/spg.even.100k"
"$GCTA_BIN" --bfile "$GENO_DIR" --make-grm --thread-num "$THREADS" --out "$PREFIX"
"$GCTA_BIN" --grm "$PREFIX" --make-bK-sparse 0.05 --thread-num "$THREADS" --out "$PREFIX.sparse"
for i in {1..10}; do
  "$GCTA_BIN" --bfile "$GENO_DIR" --fastGWA-mlm --grm-sparse "$PREFIX.sparse"     --pheno "$PHENO_DIR/pheno.$i.txt" --thread-num "$THREADS" --out "$PREFIX.pheno$i.fastgwa"
done
# GCTA appends .fastGWA: final names end in .fastgwa.fastGWA, case-sensitive.
