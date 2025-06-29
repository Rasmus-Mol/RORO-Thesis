#!/bin/sh
#BSUB -J "Deterministic_panic_test[1-8]"
#BSUB -o HPC_Output/output_%J.out
#BSUB -e HPC_Output/error_%J.err
#BSUB -q hpc
#BSUB -n 4
#BSUB -R "rusage[mem=6GB]"
#BSUB -R "span[hosts=1]"
#BSUB -W 5:00
#BSUB -u s194364@student.dtu.dk
#BSUB -N 
# end of BSUB options

# load Julia and gurobi version
#module load gurobi/1.7.0
module load julia/1.11.3
module load gurobi/12.0.1

julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'

julia Panic_Model_test.jl $LSB_JOBINDEX 
