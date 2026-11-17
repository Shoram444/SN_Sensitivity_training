#!/bin/sh

# SLURM options:
#SBATCH --partition=htc
#SBATCH --job=test
#SBATCH --licenses=sps
#SBATCH --array=0-0
#SBATCH --time=2-0
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=/sps/nemo/scratch/mpetro/Sensitivity_training/nd_%job_%A_%a.log


#################################
#################################
#################################
###### README 
###### This is an example script that shows how to run a single run of the nd optimization scrpt on CC
###### There are a few things that need to be adjusted before running the job:
###### 1. select which script to run (i.e. examples/run_nd_sensitivity_on_mathis_data.jl)
###### 2. set up the constants within the script, i.e. signal, optimization_time_limit, output_directory, etc.
###### 3. set up accordingly job parameters in SLURM
######
#################################
#################################
#################################



PROJECT_DIR=/sps/nemo/scratch/mpetro/Sensitivity_training/ # path to where Project.toml file is located (Project.toml file contains all the dependencies for the project)
SCRIPT=/sps/nemo/scratch/mpetro/Sensitivity_training/examples/run_nd_sensitivity_on_mathis_data.jl # example run script

export JULIA_NUM_THREADS=8 # selec number of threads to use, some tasks can run in parallel so that's good
export JULIA_CPU_TARGET="generic"


echo "Running script: $SCRIPT"

module load julia 
julia --project=$PROJECT_DIR -t $JULIA_NUM_THREADS --pkgimages=no $SCRIPT 
# julia --project=$PROJECT_DIR -t $JULIA_NUM_THREADS --pkgimages=no -J $PROJECT_DIR/SNSymlink/SNPackageSysimage.so $SCRIPT 

echo "Finished running script: $SCRIPT"



