#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
    echo 'Usage: sh scripts/run_hpc.sh /absolute/path/config.yaml' >&2
    exit 2
fi

# Set before importing numerical libraries. Override intentionally in your job environment.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

echo "[$(date -u +%FT%TZ)] Validating inputs"
adni-x validate --config "$1"
echo "[$(date -u +%FT%TZ)] Preparing data"
adni-x prepare --config "$1"
echo "[$(date -u +%FT%TZ)] Running analyses"
adni-x run all --config "$1"

