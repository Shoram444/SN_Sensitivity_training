#!/bin/bash
#SBATCH --job-name=sensitivity-smoke
#SBATCH --partition=htc
#SBATCH --licenses=sps
#SBATCH --time=00:10:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=smoke-%j.log

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:?Submit from the repository root using sbatch scripts/smoke_test.sh}"
PROJECT_DIR="$PWD"
if [[ ! -f "$PROJECT_DIR/Project.toml" || ! -f "$PROJECT_DIR/scripts/test_metadata.jl" ]]; then
    echo "Submit from the repository root using sbatch scripts/smoke_test.sh" >&2
    exit 1
fi

export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export JULIA_CPU_TARGET=generic
module load julia root

echo "Running metadata smoke tests in $PROJECT_DIR"
julia --project="$PROJECT_DIR" --startup-file=no --threads="$JULIA_NUM_THREADS" --pkgimages=no \
    -e 'Sys.which("root") === nothing && error("ROOT missing on compute node"); include(ARGS[1])' \
    "$PROJECT_DIR/scripts/test_metadata.jl"

echo "Running synthetic workflow smoke tests"
julia --project="$PROJECT_DIR" --startup-file=no --threads="$JULIA_NUM_THREADS" --pkgimages=no \
    "$PROJECT_DIR/scripts/test_workflow.jl"

echo "Smoke test completed successfully."
