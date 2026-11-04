#!/usr/bin/env bash
# Generic SLURM template: no guessed account, partition, memory or time.
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Template requires a confirmed SLURM allocation" >&2; exit 2; }
: "${SLURM_CPUS_PER_TASK:?Supply --cpus-per-task at submission}"
[[ $# -eq 1 && -f "$1" ]] || { echo "Pass one reviewed shell runner" >&2; exit 2; }
export THREADS="$SLURM_CPUS_PER_TASK"
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS=1
bash "$1"
