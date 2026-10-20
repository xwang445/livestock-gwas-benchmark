#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
: "${GENO_SOURCE:?Set GENO_SOURCE to the original PLINK genotype prefix}"
: "${REGENIE_GENO_SOURCE:?Set REGENIE_GENO_SOURCE}"
: "${SNP_INFO_SOURCE:?Set SNP_INFO_SOURCE}"
export OPENBLAS_NUM_THREADS=16 OMP_NUM_THREADS=16 MKL_NUM_THREADS=16
export TMPDIR="$PWD/tmp"
exec 9>status/prepare.lock
flock -n 9 || exit 1
trap 'echo "FAILED preparation at line $LINENO" > status/prepare.state' ERR
echo "RUNNING" > status/prepare.state
Rscript scripts/audit.R > logs/audit.log 2>&1
if [ ! -f status/genotype.PASS ]; then
  "${PLINK_BIN:-plink}" --bfile "$GENO_SOURCE" --keep genotype/keep.txt --keep-allele-order --make-bed --threads 8 --memory 8192 --out genotype/spg.20k > logs/genotype.log 2>&1
  "${PLINK_BIN:-plink}" --bfile genotype/spg.20k --keep-allele-order --freq --missing --threads 8 --memory 8192 --out genotype/spg.20k > logs/genotype_qc.log 2>&1
  Rscript scripts/check_genotype.R > logs/genotype_check.log 2>&1
fi
if [ ! -f status/grm.PASS ]; then
  "${MPH_BIN:-mph}" --make_grm --binary_genotype genotype/spg.20k --snp_info grm/snp_info.csv --num_threads 16 --out grm/spg.20k > logs/grm.log 2>&1
  /usr/bin/python3 scripts/check_grm.py
fi
echo "COMPLETE" > status/prepare.state
