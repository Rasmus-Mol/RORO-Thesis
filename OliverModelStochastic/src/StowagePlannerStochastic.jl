module StowagePlannerStochastic

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
# New
using Statistics
using Distributions
using HypothesisTests 
using StructTypes
using Hungarian
using Clustering

const GRB_ENV = Ref{Gurobi.Env}()
function __init__()
    const GRB_ENV[] = Gurobi.Env()
    return
end

include("representation/cargo.jl")
include("representation/deck.jl")
include("representation/slot.jl")
include("representation/vessel.jl")
include("representation/instance.jl")
include("representation/problem.jl")

include("model/base_model.jl")
#include("model/hazardous.jl")
include("model/stability.jl")

include("utils/helpers.jl")
include("solution.jl")
include("plots/solution.jl")

# New
include("model/base_stochastic_model.jl")
include("model/stability_stochastic.jl")
include("CompareSolutions.jl")
include("plots/weight_plots.jl")
include("representation/CargoScenarios.jl")
include("model/second_stage_model.jl")
include("utils/SaveData.jl")
include("representation/VarianceOfWeight.jl")
include("model/random_stowage_plan.jl")
include("representation/ScenarioReduction.jl")
include("utils/addnoise.jl")

#include("loadmaster/types.jl")
#include("loadmaster/api.jl")


# Makes these functions available outside the module
# load data and solve model function - Deterministic
export load_data, create_model, solve_model, extract_solution, plot_solution
# Stochastic functions
export create_model_stochastic, extract_stochastic_solution, plot_stochastic_solution
export create_stochastic_problem, create_stochastic_model, Monto_Carlo_sampling
export compare_solutions_print, plot_ballast_weight_diff, plot_ballast_cargo_weight_diff
export second_stage_model_slack
# Struct for stowage problem and solution
export StowageProblem, Solution
export get_solution_second_stage_stochastic, get_solution_second_stage_deterministic, second_stage_model
# Stochastic structs
export StochasticStowageProblem, SolutionStochastic
# plots
export plot_cargo_weights, plot_cargo_OG
# EVP
export expected_value_problem

# Save data to JSON files
export write_problem, write_problem_stochastic, write_solution, write_solution_stochastic, write_HPC_data
# Get data from JSON files
export get_deterministic_problem, get_stochastic_problem, get_solution_deterministic, get_solution_stochastic, get_HPC_data
export now, format, Dates, write_slack, get_slack
# Bootstrap method 
export Bootstrap_bookedweight_quantile

# Other models
export random_stowage_plan, random_cargocollection, sort_cargocollection
export create_random_stowageplan_model, create_model_stochastic_cargo_fraction
export create_random_stowageplan_model, random_cargocollection, random_stowage_plan

# Scenario reduction problem
export create_stochastic_problem_scenarioreduction, generate_simple_cargo_scenarios
export scenario_reduction_heuristic, scenario_reduction_naive
export scenario_reduced, scenario_reduction_clustering

# Noise
export add_white_noise_to_test_instance, add_biased_noise_to_test_instance

# slack
export create_model_slack_deck, create_model_stochastic_slack_deck, second_stage_model_deckslacked
export second_stage_model_deckslacked_slack

# panic model tests
export create_model_temp, calculate_vcg_slopes_temp, calculate_vcg_slopes, create_model_test
export calculate_vcg_slopes_1

# Functions to connect to LoadMaster
#export LoadMaster, getcurrentvessel, resetvessel, upload_solution, getresults, parse_result, uri, headers, struct_to_json

end # module StowagePlannerStochastic
