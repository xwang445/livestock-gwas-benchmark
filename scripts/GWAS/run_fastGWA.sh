#!/bin/bash
set -e
set -o pipefail

# =============================
# Define paths
# =============================

GENO_DIR=/path/to/genotype/spg.even.100k
PHENO_DIR=/path/to/phenotype/fastGWA_data
OUTPUT_DIR=./results
PREFIX=spg.even.100k
THREADS=16

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# =============================
# Step 1: Make GRM
# =============================

gcta_1.94.1 \
  --bfile ${GENO_DIR} \
  --make-grm \
  --thread-num "$THREADS" \
  --out ${PREFIX}

# =============================
# Step 2: Make Sparse GRM
# =============================

gcta_1.94.1 \
  --grm ${PREFIX} \
  --make-bK-sparse 0.05 \
  --thread-num $THREADS \
  --out ${PREFIX}.sparse

# =============================
# Step 3: Run fastGWA
# =============================

for pheno_file in ${PHENO_DIR}/pheno.{1..10}.txt; do

    trait=$(basename "${pheno_file}" | sed -nE 's/^pheno\.([0-9]+)\.txt$/\1/p')

    gcta_1.94.1 \
        --bfile ${GENO_DIR} \
        --fastGWA-mlm \
        --grm-sparse ${OUTPUT_DIR}/${PREFIX}.sparse \
        --pheno ${pheno_file} \
        --thread-num ${THREADS} \
        --out ${OUTPUT_DIR}/${PREFIX}.pheno${trait}.fastgwa

    echo "Finished pheno.${trait}"
    echo ""

done
