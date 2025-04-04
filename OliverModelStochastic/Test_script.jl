# Test place
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
include("src/representation/VarianceOfWeight.jl")


# load data
problem = load_data("finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous")
scenarios = 20
n_cargo_unknownweight = length(problem.cargo) # all cargo weights are unknown
stochastic_problem1 = create_stochastic_problem(problem, scenarios, n_cargo_unknownweight,[])
EVP_problem1 = expected_value_problem(stochastic_problem1)
stochastic_problem2 = create_stochastic_problem(problem, scenarios, n_cargo_unknownweight,[],Bootstrap_bookedweight_quantile)
EVP_problem2 = expected_value_problem(stochastic_problem2)

# Plot weight distribution
include("src/plots/weight_plots.jl")
#plot_weights(stochastic_problem1,1) # plot for scenario 1
plots_gen = plot_cargo_weights(stochastic_problem1,[1,2,3]) # plot for multiple scenarios
plots_boot = plot_cargo_weights(stochastic_problem2,[1,2,3])
for i in 1:length(plots_gen)
    display(plots_gen[i])
end
OG_plot = plot_cargo_OG(problem) # plot for original problem
EVP_plot = plot_cargo_OG(EVP_problem)
display(OG_plot[1])
display(OG_plot[2])
display(EVP_plot[1])
display(EVP_plot[2])
for i in 1:length(plots_boot)
    display(plots_boot[i])
end
# total weight of cargo
println("############################")
println("Total weight of cargo: ", sum([cargo.weight for cargo in problem.cargo]))
println("Average total weight of total cargo, stochastic problem - Gen: ", sum(sum([cargo.weight for cargo in stochastic_problem1.cargo.items[i]]) for i in 1:scenarios)/scenarios)
println("Average total weight of total cargo, stochastic problem - Boot: ", sum(sum([cargo.weight for cargo in stochastic_problem2.cargo.items[i]]) for i in 1:scenarios)/scenarios)
println("############################")
# Print mean weight of each cargo type
for j in 1:4
    if length([cargo.weight for cargo in filter(x -> x.cargo_type_id == j,problem.cargo)]) > 0
        println("Mean of cargo-type ", j, ": ", mean([cargo.weight for cargo in filter(x -> x.cargo_type_id == j,problem.cargo)]))
    end
end
for i in 1:scenarios
    println("############################")
    for j in 1:4
        if length([cargo.weight for cargo in filter(x -> x.cargo_type_id == j,stochastic_problem1.cargo.items[i])]) > 0
            println("Mean of cargo-type ", j, ": ", mean([cargo.weight for cargo in filter(x -> x.cargo_type_id == j,stochastic_problem1.cargo.items[i])]))
        end
    end
    println("Total weight of cargo - Gen, scenario ", i, ": ", sum([cargo.weight for cargo in stochastic_problem1.cargo.items[i]]))
end
println("############################")
for i in 1:scenarios
    println("############################")
    for j in 1:4
        if length([cargo.weight for cargo in filter(x -> x.cargo_type_id == j,stochastic_problem2.cargo.items[i])]) > 0
            println("Mean of cargo-type ", j, ": ", mean([cargo.weight for cargo in filter(x -> x.cargo_type_id == j,stochastic_problem2.cargo.items[i])]))
        end
    end
    println("Total weight of cargo - Boot, scenario ", i, ": ", sum([cargo.weight for cargo in stochastic_problem2.cargo.items[i]]))
end
println("############################")

# Model scripts
include("src/model/base_model.jl")
include("src/model/stability.jl")
include("src/model/base_stochastic_model.jl")
include("src/model/stability_stochastic.jl")
include("src/utils/helpers.jl")

# create model
model_stochastic1 = create_model_stochastic(stochastic_problem1)
model_stochastic2 = create_model_stochastic(stochastic_problem2)
model_deterministic = create_model(problem)
model_EVP1 = create_model(EVP_problem1)
model_EVP2 = create_model(EVP_problem2)

# solve model
set_silent(model_stochastic1) # removes terminal output
set_time_limit_sec(model_stochastic1, 60 * 5) # 5 minutes to solve model
#set_time_limit_sec(model_stochastic, 60 * 15) # 10 minutes to solve model
optimize!(model_stochastic1)
set_silent(model_stochastic2)
set_time_limit_sec(model_stochastic2, 60 * 5)
optimize!(model_stochastic2)
set_silent(model_deterministic) # removes terminal output
set_time_limit_sec(model_deterministic, 60 * 5) # 5 minutes to solve model
optimize!(model_deterministic)

set_silent(model_EVP1) # removes terminal output
set_time_limit_sec(model_EVP1, 60 * 5) # 5 minutes to solve model
optimize!(model_EVP1)
set_silent(model_EVP2) # removes terminal output
set_time_limit_sec(model_EVP2, 60 * 5) # 5 minutes to solve model
optimize!(model_EVP2)

# File for solution functions
include("src/solution.jl")
# Get solution
solution_stochastic1 = extract_stochastic_solution(stochastic_problem1, model_stochastic1)
solution_stochastic2 = extract_stochastic_solution(stochastic_problem2, model_stochastic2)
solution_deterministic = extract_solution(problem, model_deterministic)
solution_EVP1 = extract_solution(problem, model_EVP1)
solution_EVP2 = extract_solution(problem, model_EVP2)
println("-----------------------------------")
println("Stochastic model time, default scenarios: ", solution_stochastic1.time)
println("Stochastic model cargo loaded, default scenarios: ", solution_stochastic1.n_cargo_loaded)
println("Stochastic model status, default scenarios: ", solution_stochastic1.status)
println("Stochastic model gap, default scenarios: ", solution_stochastic1.gap)
println("Stochastic model time, Historic data scenarios: ", solution_stochastic2.time)
println("Stochastic model cargo loaded, Historic data scenarios: ", solution_stochastic2.n_cargo_loaded)
println("Stochastic model status, Historic data scenarios: ", solution_stochastic2.status)
println("Stochastic model gap, Historic data scenarios: ", solution_stochastic2.gap)
println("-----------------------------------")
println("Deterministic model time: ", solution_deterministic.time)
println("Deterministic model cargo loaded: ", solution_deterministic.n_cargo_loaded)
println("Deterministic model status: ", solution_deterministic.status)
println("-----------------------------------")
println("EVP model 1 time: ", solution_EVP1.time)
println("EVP model 1 cargo loaded: ", solution_EVP1.n_cargo_loaded)
println("EVP model 1 status: ", solution_EVP1.status)
println("EVP model 1 gap: ", solution_EVP1.gap)
println("-----------------------------------")
println("EVP model 2 time: ", solution_EVP2.time)
println("EVP model 2 cargo loaded: ", solution_EVP2.n_cargo_loaded)
println("EVP model 2 status: ", solution_EVP2.status)
println("EVP model 2 gap: ", solution_EVP2.gap)
println("-----------------------------------")

# Get solution, when knowing all cargo weights
include("src/model/second_stage_model.jl")
include("src/solution.jl")
cs_gen = solution_stochastic1.cs
cs_boot = solution_stochastic2.cs
new_model_gen = second_stage_model(cs_gen, problem)
new_model_boot = second_stage_model(cs_boot, problem)
set_silent(new_model_gen) # removes terminal output
set_silent(new_model_boot) # removes terminal output
set_time_limit_sec(new_model_gen, 60 * 5) # 5 minutes to solve model
set_time_limit_sec(new_model_boot, 60 * 5) # 5 minutes to solve model
optimize!(new_model_gen)
optimize!(new_model_boot)
fitted_solution_gen = get_solution_second_stage_stochastic(problem, new_model_gen, solution_stochastic1)
fitted_solution_boot = get_solution_second_stage_stochastic(problem, new_model_boot, solution_stochastic2)
# Do something with this information
println("-----------------------------------")
println("Second stage model cargo loaded, Gen: ", fitted_solution_gen.n_cargo_loaded)
println("Second stage model ballast weight, Gen: ", fitted_solution_gen.ballast_weight)
println("Second stage model cargo weight, Gen: ", fitted_solution_gen.cargo_weight)
println("Second stage model cargo loaded, Boot: ", fitted_solution_boot.n_cargo_loaded)
println("Second stage model ballast weight, Boot: ", fitted_solution_boot.ballast_weight)
println("Second stage model cargo weight, Boot: ", fitted_solution_boot.cargo_weight)
println("-----------------------------------")
# only works ids is numerated from 1
#=weight_diff = 0
for i in 1:length(problem.cargo)
    if sum(cs_gen,dims = 2)[i] > sum(cs_monte,dims=2)[i]
        println("Cargo ", i, " has different placements")
        weight_diff += problem.cargo[i].weight
    elseif sum(cs_gen,dims=2)[i] < sum(cs_monte,dims=2)[i]
        println("Cargo ", i, " has different placements")
        weight_diff -= problem.cargo[i].weight
    end
end
println("Total weight difference: ", weight_diff)
println("Check: ", fitted_solution_gen.cargo_weight - fitted_solution_monte.cargo_weight)
println("-----------------------------------")
=#

# Plot solution
include("src/plots/solution.jl")
plot_solution(solution_deterministic)
plot_solution_stochastic(solution_stochastic1,1)
plot_solution(solution_EVP1)
plot_solution_stochastic(solution_stochastic2,1)
plot_solution(solution_EVP2)

# Compare / print
include("src/CompareSolutions.jl")
diff_ballast, diff_cargo_weight, diff_placed, diff_index, same_index = compare_solutions_print(solution_stochastic2, solution_deterministic,1)
println("Number of different placements: ", diff_placed)
plot_ballast_weight_diff(diff_ballast,scenarios,n_cargo_unknownweight)
