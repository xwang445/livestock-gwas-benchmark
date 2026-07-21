#!/bin/bash

# =============================
# Define paths
# =============================

GENO_DIR=/path/to/genotype/spg.even.100k
PHENO_FILE=/path/to/phenotype/pheno_combined.txt
OUTPUT_DIR=./results
THREADS=16
BSIZE=1000

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# =============================
# Run REGENIE (without LOCO) for each phenotype
# =============================

for i in {1..10}
do
  regenie \
    --step 1 \
    --bed ${GENO_DIR} \
    --phenoFile ${PHENO_FILE} \
    --phenoCol pheno$i \
    --print-prs \
    --qt \
    --bsize ${BSIZE} \
    --threads ${THREADS} \
    --out regenie_step1_pheno$i

  regenie \
    --step 2 \
    --bed ${GENO_DIR} \
    --phenoFile ${PHENO_FILE} \
    --phenoCol pheno$i \
    --qt \
    --use-prs \
    --pred regenie_step1_pheno$i"_prs.list" \
    --bsize ${BSIZE} \
    --threads ${THREADS} \
    --out regenie_pheno_wo_loco$i
done
