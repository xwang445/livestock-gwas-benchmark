#!/usr/bin/env bash
# Run L2 (LOCO GRM + eigendecomposition + REGENIE eta + per-locus adjusted LD) for all chromosomes.
# nohup bash loco/scripts/run_L2_all.sh > loco/logs/run_L2_all.log 2>&1 < /dev/null &
set -euo pipefail
L="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$L"
NPAR=${1:-6}; THREADS=${2:-12}
echo "$(date -Is) start NPAR=$NPAR THREADS=$THREADS"
# largest chromosomes first
printf "%s\n" 1 13 2 6 14 8 9 4 15 3 7 5 11 16 10 12 18 17 | \
  xargs -P "$NPAR" -I{} bash -c "OPENBLAS_NUM_THREADS=$THREADS OMP_NUM_THREADS=$THREADS python3 scripts/L2_loco_ld.py {} > logs/L2_chr{}.log 2>&1"
echo "$(date -Is) finished: DONE=$(ls status/L2_chr*.DONE 2>/dev/null | wc -l)/18"
