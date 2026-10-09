#!/bin/bash
#SBATCH --job-name=adni-x
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=adni-x-%j.log

# Submit from the repository root. Resources are initial pilot estimates only.
# Usage: sbatch scripts/qmul_slurm.sh validate --config config.local.yaml
#        sbatch scripts/qmul_slurm.sh run all --config config.pilot.yaml --include-sensitivity
#        sbatch scripts/qmul_slurm.sh pipeline /absolute/path/config.local.yaml
set -eo pipefail

if [ -z "${SLURM_JOB_ID:-}" ]; then
    echo 'Submit this script with sbatch from the repository root.' >&2
    exit 2
fi
if [ "$#" -eq 0 ]; then
    echo 'Supply adni-x arguments, or pipeline CONFIG.' >&2
    exit 2
fi
cd "${SLURM_SUBMIT_DIR:?Missing Slurm submission directory}"
if [ ! -f pyproject.toml ] || [ ! -f scripts/run_hpc.sh ]; then
    echo 'Submit from the ADNI-X repository root.' >&2
    exit 2
fi
: "${ADNI_X_PYTHON_MODULE:?Export the exact Python module used to create the virtual environment}"
module load "$ADNI_X_PYTHON_MODULE"
source "${ADNI_X_VENV:-$PWD/.venv}/bin/activate"
set -u

# One process; prevent inherited BLAS settings from oversubscribing the job.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export NUMEXPR_NUM_THREADS="$OMP_NUM_THREADS"
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1

echo "Started $(date -u +%FT%TZ) on $(hostname); job $SLURM_JOB_ID"
python --version
python -m pip check
adni-x --version

if [ "$1" = pipeline ]; then
    shift
    bash scripts/run_hpc.sh "$@"
else
    adni-x "$@"
fi
echo "Completed $(date -u +%FT%TZ)"
