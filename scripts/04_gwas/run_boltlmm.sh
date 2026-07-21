#!/bin/bash

# =============================
# Define paths
# =============================

GENO_DIR=/path/to/genotype/spg.even.100k
PHENO_FILE=/path/to/phenotype/pheno_combined.txt
OUTPUT_DIR=./results
THREADS=16

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# =============================
# Run BOLT-LMM for each phenotype
# =============================

for i in {1..10}
do
  bolt \
  --bfile ${GENO_DIR} \
  --phenoFile ${PHENO_FILE} \
  --phenoCol pheno$i \
  --lmm \
  --LDscoresUseChip \
  --numThreads ${THREADS} \
  --statsFile bolt_pheno$i.stats
done
