push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
#include("src/representation/ScenarioReduction.jl")
using .StowagePlannerStochastic
using JuMP
HPC_folder = "HPC_test"
if !isdir("Results/"*HPC_folder)
    mkdir("Results/"*HPC_folder)
end

temp = rand(1:10,5,5)
nested = [[temp[i, j] for j in 1:size(temp, 2)] for i in 1:size(temp, 1)]
open(joinpath("Results",HPC_folder,"test"*".json"), "w") do file
    JSON.print(file, nested, 4)
end
