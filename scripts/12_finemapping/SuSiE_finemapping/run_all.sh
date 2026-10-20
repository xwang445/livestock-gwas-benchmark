#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 scripts/01_ordinary_ld.py > logs/ordinary_ld.log 2>&1
Rscript scripts/02_susie_rss.R 20 > logs/susie_rss.log 2>&1
python3 scripts/03_evaluate_and_plot.py > logs/evaluation.log 2>&1
touch status/RUN.COMPLETE
