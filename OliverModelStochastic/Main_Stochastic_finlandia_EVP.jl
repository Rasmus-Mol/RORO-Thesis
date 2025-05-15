# Script to run on HPC - Runs EVP problems
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
using .StowagePlannerStochastic
using JuMP