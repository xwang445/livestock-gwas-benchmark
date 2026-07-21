#!/bin/bash
set -e
set -o pipefail

########################################
# Configuration
########################################

THREADS=16
export OMP_NUM_THREADS=${THREADS}
GENO_DIR=/path/to/genotype/spg.even.100k
PHENO_FILE=/path/to/phenotype/pheno.slemm.csv
SNP_INFO_FILE=/path/to/genotype/snp_info.csv
OUTPUT_DIR=./results

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# Set local aliases for compatibility with script
BFILE=${GENO_DIR}
PFILE=${GENO_DIR}
PHENO=${PHENO_FILE}
SNP_INFO=${SNP_INFO_FILE}

########################################
# Loop over traits
########################################

for i in $(seq 1 10); do

  TRAIT=pheno${i}
  PREFIX=spg.even.100k.${TRAIT}

  echo "======================================"
  echo "Running GWAS for ${TRAIT}"
  echo "======================================"

  ########################################
  # Step 1 — Fit null LMM
  ########################################

  echo "Step 1: Fitting null LMM..."

  slemm \
    --lmm \
    --phenotype_file "${PHENO}" \
    --trait ${TRAIT} \
    --bfile "${BFILE}" \
    --snp_info_file "${SNP_INFO}" \
    --out ${PREFIX}

  ########################################
  # Step 2 — Compute GRAMMAR-Gamma
  ########################################

  echo "Step 2: Computing gamma..."

  OMP_NUM_THREADS=1 \
  slemm_gamma \
    --pfile "${PFILE}" \
    --slemm ${PREFIX} \
    --out ${PREFIX}.gamma.txt

  ########################################
  # Step 3 — GWAS per chromosome
  ########################################

  echo "Step 3: Running per-chromosome GWAS..."

  for chr in $(seq 1 18); do
    slemm_gwa \
      --pfile "${PFILE}" \
      --slemm ${PREFIX} \
      --chr ${chr} \
      --out ${PREFIX}.chr${chr}.gwas.txt
  done

  ########################################
  # Step 4 — Merge results
  ########################################

  echo "Merging chromosome results..."

  cp ${PREFIX}.chr1.gwas.txt ${PREFIX}.allchr.gwas.txt

  for chr in $(seq 2 18); do
    tail -n +2 ${PREFIX}.chr${chr}.gwas.txt >> ${PREFIX}.allchr.gwas.txt
    rm ${PREFIX}.chr${chr}.gwas.txt
  done

  echo "Finished ${TRAIT}"
  echo ""

done

echo "All traits completed successfully."
