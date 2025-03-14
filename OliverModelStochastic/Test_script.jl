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

include("src/representation/cargo.jl")
include("src/representation/deck.jl")
include("src/representation/slot.jl")
include("src/representation/vessel.jl")
include("src/representation/instance.jl")
include("src/representation/problem.jl")
# New: Stochastic 
include("src/representation/CargoScenarios.jl")


# load data
scenarios = 10
n_cargo_unknownweight = 50
problem = load_data("finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous")
stochastic_problem = create_stochastic_problem(problem, scenarios, n_cargo_unknownweight)

# Model scripts
include("src/model/base_model.jl")
include("src/model/stability.jl")
include("src/model/base_stochastic_model.jl")
include("src/model/stability_stochastic.jl")
include("src/utils/helpers.jl")

# create model
model_stochastic = create_model_stochastic(stochastic_problem)
model_deterministic = create_model(problem)

# solve model
set_silent(model_stochastic) # removes terminal output
set_time_limit_sec(model_stochastic, 60 * 5) # 5 minutes to solve model
#set_time_limit_sec(model_stochastic, 60 * 10) # 10 minutes to solve model
optimize!(model_stochastic)
set_silent(model_deterministic) # removes terminal output
set_time_limit_sec(model_deterministic, 60 * 5) # 5 minutes to solve model
optimize!(model_deterministic)

# File for solution functions
include("src/solution.jl")

# Get solution
solution_stochastic = extract_stochastic_solution(stochastic_problem, model_stochastic)
solution_deterministic = extract_solution(problem, model_deterministic)

println("Stochastic model time: ", solution_stochastic.time)
println("Stochastic model cargo loaded: ", solution_stochastic.n_cargo_loaded)
println("Deterministic model time: ", solution_deterministic.time)
println("Deterministic model cargo loaded: ", solution_deterministic.n_cargo_loaded)


# Plot solution
include("src/plots/solution.jl")
plot_solution(solution_deterministic)
plot_solution_stochastic(solution_stochastic,1)

# Compare / print
include("src/CompareSolutions.jl")
diff_ballast = compare_solutions_print(solution_stochastic, solution_deterministic,1)
plot_ballast_weight_diff(diff_ballast,scenarios,n_cargo_unknownweight)


# Test different amount of scenarios and different number of unknown cargo cargo_weights
# Takes a long time to run
#=
scenarios = [5,10,20,50]
n_unknown_weights = [5,20,50,100]
model_deterministic = create_model(problem)
set_silent(model_deterministic) # removes terminal output
set_time_limit_sec(model_deterministic, 60 * 5) # 5 minutes to solve model
optimize!(model_deterministic)
solution_deterministic = extract_solution(problem, model_deterministic)
s_solutions = []
s_problems = []
for i in 1:length(scenarios)
    for j in 1:length(n_unknown_weights)
        println("Scenarios: ", scenarios[i], ". Number of unknown weights: ", n_unknown_weights[j])
        stochastic_problem = create_stochastic_problem(problem, scenarios[i], n_unknown_weights[j])
        s_problems = push!(s_problems,stochastic_problem)
        model_stochastic = create_model_stochastic(stochastic_problem)
        set_silent(model_stochastic) # removes terminal output
        set_time_limit_sec(model_stochastic, 60 * 5) # 5 minutes to solve model
        optimize!(model_stochastic)
        solution_stochastic = extract_stochastic_solution(stochastic_problem, model_stochastic)
        push!(s_solutions, solution_stochastic)
        diff_ballast,diff_cargo = compare_solutions_print(solution_stochastic, solution_deterministic,1)
        plot_ballast_weight_diff(diff_ballast,scenarios[i],n_unknown_weights[j])
    end
end
=#
diff_ballast, diff_cargo = compare_solutions_print(s_solutions[1], solution_deterministic,1)
plot_ballast_weight_diff(diff_ballast,scenarios[1],n_unknown_weights[1])
plot_ballast_cargo_weight_diff(diff_ballast, diff_cargo,scenarios[1],n_unknown_weights[1])

# Plot weight distribution
include("src/plots/weight_plots.jl")
plot_weights(stochastic_problem,1) # plot for scenario 1
plots = plot_cargo_weights(stochastic_problem,[1,2,3]) # plot for multiple scenarios
for i in 1:length(plots)
    display(plots[i])
end
OG_plot = plot_cargo_OG(problem) # plot for original problem







