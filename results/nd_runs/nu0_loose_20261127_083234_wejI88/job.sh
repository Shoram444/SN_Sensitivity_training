#!/usr/bin/env bash
set -euo pipefail
export SN_PROJECT_DIR=/Users/maros.petro/Work/Phd/SN_Sensitivity_training
export SN_OUTPUT_PATH=/Users/maros.petro/Work/Phd/SN_Sensitivity_training/results/nd_runs/nu0_loose_20261127_083234_wejI88/results.csv
export SN_SAVE_AUXILIARY_DATA=true
export JULIA_NUM_THREADS=6
export JULIA_CPU_TARGET=generic
cd /Users/maros.petro/Work/Phd/SN_Sensitivity_training
# module load julia
julia --project=/Users/maros.petro/Work/Phd/SN_Sensitivity_training -t 6 --pkgimages=no -J /Users/maros.petro/Work/Phd/SN_Sensitivity_training/SNSymlink/SNPackageSysimage.so /Users/maros.petro/Work/Phd/SN_Sensitivity_training/results/nd_runs/nu0_loose_20261127_083234_wejI88/analysis.jl
exec julia --project=/Users/maros.petro/Work/Phd/SN_Sensitivity_training -t 6 --pkgimages=no -J /Users/maros.petro/Work/Phd/SN_Sensitivity_training/SNSymlink/SNPackageSysimage.so /Users/maros.petro/Work/Phd/SN_Sensitivity_training/scripts/compare_data_simulation.jl --results=/Users/maros.petro/Work/Phd/SN_Sensitivity_training/results/nd_runs/nu0_loose_20261127_083234_wejI88/results.csv
