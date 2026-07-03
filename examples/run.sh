SBATCH time=48h

module load julia 
julia --project="sps_mpetro/Sensitivity_training/" sps_mpetro/Sensitivity_training/examples/ex4_freqND.jl 