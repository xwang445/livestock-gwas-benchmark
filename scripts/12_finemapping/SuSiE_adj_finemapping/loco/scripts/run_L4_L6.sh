#!/usr/bin/env bash
# After L2 (all 18 chromosomes): SuSiE-adj on every locus (L4), evaluation (L5), figures (L6). Restartable.
# nohup bash loco/scripts/run_L4_L6.sh > loco/logs/run_L4_L6.log 2>&1 < /dev/null &
set -eu
L="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$L"
log() { echo "$(date -Is) $*"; }
n=$(ls status/L2_chr*.DONE | wc -l); log "L2 DONE chromosomes: $n/18"
[ "$n" -eq 18 ] || { log "L2 incomplete - stop"; exit 1; }
n_ld=$(find ld -name summary.tsv | wc -l); log "LD summaries: $n_ld (manifest 13301)"
log "L4 SuSiE-adj (ERV=TRUE primary)"; OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript scripts/L4_susie_adj_loco.R 40 > logs/L4_susie_adj_loco.log 2>&1
log "L5 evaluation";                   OPENBLAS_NUM_THREADS=8 python3 scripts/L5_evaluate_loco.py > logs/L5_evaluate_loco.log 2>&1
log "L6 figures";                      python3 scripts/L6_figures_loco.py > logs/L6_figures_loco.log 2>&1
log "COMPLETE"; touch status/run_L4_L6.COMPLETE
