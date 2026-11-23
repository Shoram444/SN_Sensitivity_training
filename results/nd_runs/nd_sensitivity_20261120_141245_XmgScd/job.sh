#!/usr/bin/env bash
set -euo pipefail
export SN_PROJECT_DIR=/Users/maros.petro/Work/Phd/SN_Sensitivity_training
export SN_OUTPUT_PATH=/Users/maros.petro/Work/Phd/SN_Sensitivity_training/results/nd_runs/nd_sensitivity_20261120_141245_XmgScd/results.csv
export JULIA_NUM_THREADS=6
export JULIA_CPU_TARGET=generic
cd /Users/maros.petro/Work/Phd/SN_Sensitivity_training
# module load julia
exec julia --project=/Users/maros.petro/Work/Phd/SN_Sensitivity_training -t 6 --pkgimages=no -J /Users/maros.petro/Work/Phd/SN_Sensitivity_training/SNSymlink/SNPackageSysimage.so /Users/maros.petro/Work/Phd/SN_Sensitivity_training/results/nd_runs/nd_sensitivity_20261120_141245_XmgScd/analysis.jl
