#!/bin/sh
#BSUB -J "Stochastic_Model_Tests"
#BSUB -o HPC_Output/output_%J.out
#BSUB -q hpc
#BSUB -n 4
#BSUB -R "rusage[mem=2GB]"
#BSUB -R "span[hosts=1]"
#BSUB -W 10:00
#BSUB -u s194364@student.dtu.dk
#BSUB -N 
# end of BSUB options

# load Julia version

julia -e using Pkg; Pkg.activate(".")
julia -e using Pkg; pkg.instantiate()

julia Main_Stochastic.jl