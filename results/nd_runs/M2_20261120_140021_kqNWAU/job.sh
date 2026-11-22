#!/usr/bin/env bash
set -euo pipefail
export SN_PROJECT_DIR=/sps/nemo/scratch/mpetro/Sensitivity_training
export SN_OUTPUT_PATH=/sps/nemo/scratch/mpetro/Sensitivity_training/results/nd_runs/M2_20261120_140021_kqNWAU/results.csv
export JULIA_NUM_THREADS=6
export JULIA_CPU_TARGET=generic
cd /sps/nemo/scratch/mpetro/Sensitivity_training
module load julia
exec julia --project=/sps/nemo/scratch/mpetro/Sensitivity_training -t 6 --pkgimages=no -J /sps/nemo/scratch/mpetro/Sensitivity_training/SNSymlink/SNPackageSysimage.so /sps/nemo/scratch/mpetro/Sensitivity_training/results/nd_runs/M2_20261120_140021_kqNWAU/analysis.jl
