#!/bin/sh
#BSUB -J "Stochastic_model_slackdeck_instance_7"
#BSUB -o HPC_Output/output_%J.out
#BSUB -e HPC_Output/error_%J.err
#BSUB -q hpc
#BSUB -n 4
#BSUB -R "rusage[mem=5GB]"
#BSUB -R "span[hosts=1]"
#BSUB -W 50:00
#BSUB -u s194364@student.dtu.dk
#BSUB -N 
# end of BSUB options

# load Julia and gurobi version
#module load gurobi/1.7.0
module load julia/1.11.3
module load gurobi/12.0.1

julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'

julia Main_Stochastic_finlandia_slacked_deck_limit.jl
