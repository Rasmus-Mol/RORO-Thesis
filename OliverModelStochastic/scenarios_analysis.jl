# Script to plot data from problem instances and historic data 
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
using StatsPlots

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
include("src/representation/VarianceOfWeight.jl")
# Test instances
include("src/utils/test_instances.jl")


# Scenario analysis
plot_folder = "Plots/Data/Scenarios/"
# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end
# Load data - change this to be correct
HPC_folder = "Finlandia_no_cars_light_6014_05_20"
repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
sc = length(scenarios)
n = length(n_unknown)
test_problem = Finlandia_test[5]
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", test_problem, "hazardous"
Deterministic_problem = load_data(problemname1,problemname2,problemname3)
EVP_problem_gen = Array{Any}(undef, repetitions, sc,n)
Stochastic_problem_gen = Array{Any}(undef, repetitions, sc,n)
EVP_problem_boot = Array{Any}(undef, repetitions, sc,n)
Stochastic_problem_boot = Array{Any}(undef, repetitions, sc,n)

for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # Load EVP
            foldername = "EVP_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Problem"
            EVP_problem_gen[i,j,k] = get_deterministic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
            foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Problem"
            EVP_problem_boot[i,j,k] = get_deterministic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
            # Load Stochastic
            foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem_gen[i,j,k] = get_stochastic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
            foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem_boot[i,j,k] = get_stochastic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
        end
    end
end
 # plot for original problem
OG_plot = plot_cargo_OG(Deterministic_problem)
display(OG_plot[1])
display(OG_plot[2])
# Save plots
savefig(OG_plot[1],plot_folder*"Cargo_distribution_Determinstic_1_"*test_problem*".png")
savefig(OG_plot[2],plot_folder*"Cargo_distribution_Determinstic_2"*test_problem*".png")
# uniform random sampling
plots_gen = plot_cargo_weights(Stochastic_problem_gen[1,1,1],[i for i =1:scenarios[1]])
xplots = 2
yplots = 2
plot(plots_gen[1],plots_gen[2],plots_gen[3],plots_gen[4], layout=(xplots,yplots))
total_weight_gen = []
weight_diff = []
det_weight = Deterministic_problem.cargo.total_weight
push!(total_weight_gen, Stochastic_problem_gen[1,1,1].cargo.items.total_weight)
push!(weight_diff,total_weight_gen[1].-det_weight)
p = plot([i for i=1:scenarios[1]], total_weight_gen[1]-ones(scenarios[1])*det_weight, label = "No. scenarios: $(scenarios[1])", xlabel = "Scenario", ylabel = "Weight (t.)", title = "Weight difference")
for i in 2:sc
    push!(total_weight_gen, Stochastic_problem_gen[1,i,1].cargo.items.total_weight)
    push!(weight_diff,total_weight_gen[i].-det_weight)
    p = plot!([j for j=1:scenarios[i]], total_weight_gen[i]-ones(scenarios[i])*det_weight, label = "No. scenarios: $(scenarios[i])")
end
display(p)
savefig(p,plot_folder*"Weight_difference_uniform_random_sampling_nocars_light60.png")
# generate scenarios for other test - should be loaded but is not made yet
gen_problem_nocars_heavy_60 = Array{Any}(undef, sc)
nocars_heavy_60 = Finlandia_test[7]
problemname1, problemname2, problemname3 = "finlandia", nocars_heavy_60, "hazardous"
Det_pro_nocars_heavy_60 = load_data(problemname1,problemname2,problemname3)
det_weight_nocars_heavy_60 = Det_pro_nocars_heavy_60.cargo.total_weight
total_weight_gen_nocars_heavy_60 = []
for i in 1:sc
    gen_problem_nocars_heavy_60[i] = create_stochastic_problem(Det_pro_nocars_heavy_60, scenarios[i],length(Det_pro_nocars_heavy_60.cargo) , []) 
    push!(total_weight_gen_nocars_heavy_60, gen_problem_nocars_heavy_60[i].cargo.items.total_weight)
end
p = plot([i for i=1:scenarios[1]], total_weight_gen_nocars_heavy_60[1]-ones(scenarios[1])*det_weight_nocars_heavy_60, label = "No. scenarios: $(scenarios[1])", xlabel = "Scenario", ylabel = "Weight (t.)", title = "Weight difference")
for i in 2:sc
    p = plot!([j for j=1:scenarios[i]], total_weight_gen_nocars_heavy_60[i]-ones(scenarios[i])*det_weight_nocars_heavy_60, label = "No. scenarios: $(scenarios[i])")
end
display(p)
savefig(p,plot_folder*"Weight_difference_uniform_random_sampling_nocars_heavy60.png")

item_id = 1
one_item_heavy = []
one_item_light = []
for i in 1:scenarios[sc]
    push!(one_item_heavy, gen_problem_nocars_heavy_60[sc].cargo.items[i].items[item_id].weight)
    push!(one_item_light, Stochastic_problem_gen[1,sc,1].cargo.items[i].items[item_id].weight)
end
p = scatter([1:scenarios[sc]], one_item_heavy, label = "Scenario weight", xlabel = "Scenario", ylabel = "Weight (t.)", title = "Weight of cargo id $(item_id)")
p = plot!([1:scenarios[sc]], ones(scenarios[sc])*Det_pro_nocars_heavy_60.cargo.items[item_id].weight, label = "Actual weight")
savefig(p,plot_folder*"Weight_uniform_random_sampling_nocars_heavy60_item_$(item_id).png")
p = scatter([1:scenarios[sc]], one_item_light, label = "Scenario weight", xlabel = "Scenario", ylabel = "Weight (t.)", title = "Weight of cargo id $(item_id)")
p = plot!([1:scenarios[sc]], ones(scenarios[sc])*Deterministic_problem.cargo.items[item_id].weight, label = "Actual weight")
savefig(p,plot_folder*"Weight_uniform_random_sampling_nocars_light60_item_$(item_id).png")


Det_pro_nocars_heavy_60.cargo.items