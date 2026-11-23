#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --job=data
#SBATCH --licenses=sps
#SBATCH --array=0-0
#SBATCH --time=0-3
#SBATCH --mem=32G
#SBATCH --cpus-per-task=6
#SBATCH --output=/sps/nemo/scratch/mpetro/Sensitivity_training/nd_%job_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=maros.petro@cvut.cz


PROJECT_DIR=/sps/nemo/scratch/mpetro/Sensitivity_training/ # path to where Project.toml file is located (Project.toml file contains all the dependencies for the project)
SCRIPT=/sps/nemo/scratch/mpetro/Sensitivity_training/scripts/compare_data_simulation.jl # example run script
RESULTS_DIR=/sps/nemo/scratch/mpetro/Sensitivity_training/results/nd_runs/M2_basic_20261123_100524_wHkw9Y

export JULIA_NUM_THREADS=6 # selec number of threads to use, some tasks can run in parallel so that's good
export JULIA_CPU_TARGET="generic"


echo "Running script: $SCRIPT"

module load julia 


cd $PROJECT_DIR

S=/sps/nemo/scratch/granjon/double_beta/cimrman_curve/new_threshold_make_root_without_ptd_simu 

{ [ -e data/mathis_files ] || ln -s $S/ data/mathis_files; }
{ [ -e data/mathis_files_reduced ] || ln -s $S/reduces_files data/mathis_files_reduced; } 

julia --project=. -t $JULIA_NUM_THREADS --pkgimages=no -J SNSymlink/SNPackageSysimage.so $SCRIPT --results=$RESULTS_DIR/results.csv 2>&1 | tee $RESULTS_DIR/compare_data.log