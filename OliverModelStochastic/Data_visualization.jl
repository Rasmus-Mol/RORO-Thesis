# Script to plot data
push!(LOAD_PATH, pwd())
using Base: @kwdef
import Base: length, getindex, lastindex, keys, eachindex, iterate, eltype, firstindex

using Random
using StatsBase
using StructArrays
using JSON3, JSON
using DataFrames
using CSV, XLSX
using Dates
using Interpolations
using UnPack
using JuMP, HiGHS, GLPK, Gurobi
using Plots
using HTTP
using URIs
# News
using Statistics
using HypothesisTests
using Distributions
using StructTypes

include("src/representation/cargo.jl")
include("src/representation/deck.jl")
include("src/representation/slot.jl")
include("src/representation/vessel.jl")
include("src/representation/instance.jl")
include("src/representation/problem.jl")
# New: Stochastic 
include("src/representation/CargoScenarios.jl")
# Plot weight distribution
include("src/plots/weight_plots.jl")
# Solution 
include("src/solution.jl")
include("src/CompareSolutions.jl")
include("src/plots/solution.jl")
# Load data Script
include("src/utils/SaveData.jl")


# load data from problems
HPC_folder = "Finlandia_27_03_13_34_52"
repetitions, scenarios, n_unknown, time_limit = get_HPC_data(HPC_folder)
repetitions = 1 # didn't finish running
sc = length(scenarios)
n = length(n_unknown)
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous"
Deterministic_problem = load_data(problemname1,problemname2,problemname3)
EVP_problem = Array{Any}(undef, repetitions, sc,n)
Stochastic_problem = Array{Any}(undef, repetitions, sc,n)

for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # Load EVP
            foldername = "EVP_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Problem"
            EVP_problem[i,j,k] = get_deterministic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
            # Load Stochastic
            foldername = "Stochastic_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem[i,j,k] = get_stochastic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
        end
    end
end

# Plot Cargo weights for some scenarios for a problem 
plots_gen = plot_cargo_weights(Stochastic_problem[1,1,5],[1,2,3]) # plot for multiple scenarios
for i in 1:length(plots_gen)
    display(plots_gen[i])
end
OG_plot = plot_cargo_OG(Deterministic_problem) # plot for original problem
EVP_plot = plot_cargo_OG(EVP_problem[1,1,1],false) # plot for EVP
display(OG_plot[1])
display(OG_plot[2])
display(EVP_plot[1])
display(EVP_plot[2])
# total weight of cargo
println("Total weight of cargo: ", sum([cargo.weight for cargo in Deterministic_problem.cargo]))
println("Average total weight of cargo, stochastic problem - Gen: ", sum(sum([cargo.weight for cargo in Stochastic_problem[1,1,1].cargo.items[i]]) for i in 1:scenarios[1])/scenarios[1])
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            println("####################")
            for l in 1:scenarios[i]
                println("Total weight of cargo - Gen, repetition: ", i, ", number of scenarios: ", scenarios[i],
                ", number of unknown weights: ", n_unknown[k], ": ",
                sum([cargo.weight for cargo in Stochastic_problem[i,j,k].cargo.items[l]]))
            end
        end
    end
end

# Plot for total weight in each scenario 
# Change this if want problems with a different number of scenarios
n_sc = 1
total_weights_sto = Array{Any}(undef, scenarios[n_sc],n)
for i in 1:n
    sto_pro = Stochastic_problem[1,n_sc,i]
    for j in 1:scenarios[n_sc]
        total_weights_sto[j,i] = sum([cargo.weight for cargo in sto_pro.cargo.items[j]])
    end
end
total_weight_det = ones(sto_pro.scenarios)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
p = plot(total_weights_sto[:,1],xlabel = "Scenario", ylabel = "Total cargo weight (t)", 
label = "Sto_pro, n: $(n_unknown[1])",
title = "Total weight for the stochastic problem,\nscenarios: $(scenarios[n_sc])")
for i in 2:n
    plot!(total_weights_sto[:,i], label = "Sto_pro, n: $(n_unknown[i])")
end
plot!(total_weight_det, label = "Det_pro")
# NB! 
# Clearly shows that our generation of weight is generally 
# higher than the deterministic problem
display(p)

# Display EVP data
n_sc = 1
total_weights_EVP = []
for i in 1:n
    push!(total_weights_EVP, sum([cargo.weight for cargo in EVP_problem[1,n_sc,i].cargo]))
end
total_weight_det = ones(n)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
p = plot(n_unknown,total_weights_EVP,xlabel = "n_unknown", ylabel = "Total cargo weight (t)", label = "EVP", 
title = "Total weight for the EVP, scenarios: $(scenarios[n_sc]), ")
plot!(n_unknown,total_weight_det, label = "Det_pro")
# Clearly shows that our generation of weight is generally 
# higher than the deterministic problem
display(p)

# Other methods for generating scenarios...


