#!/usr/bin/env bash
set -euo pipefail

# Simulate 10 polygenic phenotypes for the GWAS benchmarking study.
#
# Assumptions:
#   1. The PLINK genotype prefix contains the pig genotype data used in the study.
#   2. The SNP information file contains the fixed set of 10,000 SNPs sampled
#      from the predefined polygenic-effect windows.
#   3. MPH uses the GRM listed in GRM_LIST to simulate independent phenotypes
#      at the requested narrow-sense heritability.
#
# Usage:
#   ./simulate_polygenic_phenotypes.sh GENOTYPE_PREFIX SNP_INFO_FILE [OUT_DIR]
#
# Optional environment variables:
#   MPH_BIN, NUM_THREADS, NUM_PHENOTYPES, HERITABILITY

if (( $# < 2 || $# > 3 )); then
    echo "Usage: $0 GENOTYPE_PREFIX SNP_INFO_FILE [OUT_DIR]" >&2
    exit 2
fi

GENOTYPE_PREFIX="$1"
SNP_INFO_FILE="$2"
OUT_DIR="${3:-./simulation_output}"

MPH_BIN="${MPH_BIN:-mph}"
NUM_THREADS="${NUM_THREADS:-20}"
NUM_PHENOTYPES="${NUM_PHENOTYPES:-10}"
HERITABILITY="${HERITABILITY:-0.5}"

GRM_PREFIX="${OUT_DIR}/MPH_grm_chr1_18"
GRM_LIST="${OUT_DIR}/grmlist.txt"
PHENO_PREFIX="${OUT_DIR}/pheno"

mkdir -p "${OUT_DIR}"

# Basic input checks
for file in \
    "${GENOTYPE_PREFIX}.bed" \
    "${GENOTYPE_PREFIX}.bim" \
    "${GENOTYPE_PREFIX}.fam" \
    "${SNP_INFO_FILE}"
do
    if [[ ! -f "${file}" ]]; then
        echo "ERROR: Required input file not found: ${file}" >&2
        exit 1
    fi
done

if ! command -v "${MPH_BIN}" >/dev/null 2>&1; then
    echo "ERROR: MPH executable not found: ${MPH_BIN}" >&2
    exit 1
fi

if [[ ! "${NUM_THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: NUM_THREADS must be a positive integer." >&2
    exit 2
fi

if [[ ! "${NUM_PHENOTYPES}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: NUM_PHENOTYPES must be a positive integer." >&2
    exit 2
fi

echo "[1/3] Building the GRM from the fixed 10,000-SNP set..."

"${MPH_BIN}" \
    --binary_genotype_file "${GENOTYPE_PREFIX}" \
    --snp_info_file "${SNP_INFO_FILE}" \
    --make_grm \
    --num_threads "${NUM_THREADS}" \
    --out "${GRM_PREFIX}"

echo "[2/3] Writing the GRM list..."

# MPH expects a text file listing the GRM prefix(es), one per line.
printf '%s\n' "${GRM_PREFIX}" > "${GRM_LIST}"

echo "[3/3] Simulating ${NUM_PHENOTYPES} phenotypes at h2=${HERITABILITY}..."

"${MPH_BIN}" \
    --simulate \
    --num_threads "${NUM_THREADS}" \
    --num_phenotypes "${NUM_PHENOTYPES}" \
    --grm_list "${GRM_LIST}" \
    --heritability "${HERITABILITY}" \
    --output "${PHENO_PREFIX}"

echo "Simulation completed."
echo "GRM prefix: ${GRM_PREFIX}"
echo "Phenotype output prefix: ${PHENO_PREFIX}"
