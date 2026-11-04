#!/usr/bin/env bash
set -euo pipefail
# Recovered from ld_10Mb.log, PLINK v1.9.0-b.8 (22 Oct 2024).
# Full 31,301-animal resource; NOT the 20,000-animal GWAS subset.
: "${GENO_DIR:?Set the authorized full-resource PLINK prefix}"
: "${OUTPUT_DIR:?Set a NEW output directory}"
: "${THREADS:?Set an approved thread count}"
: "${MEMORY_MB:?Set an approved PLINK memory limit in MB}"
[[ ! -e "$OUTPUT_DIR" ]] || { echo 'Output exists' >&2; exit 2; }
mkdir -p "$OUTPUT_DIR"
"${PLINK_BIN:-plink}" --bfile "$GENO_DIR" --r2 gz --ld-window 999999 \
  --ld-window-kb 10000 --ld-window-r2 0 --seed 1773948362 \
  --threads "$THREADS" --memory "$MEMORY_MB" --out "$OUTPUT_DIR/ld_10Mb"
