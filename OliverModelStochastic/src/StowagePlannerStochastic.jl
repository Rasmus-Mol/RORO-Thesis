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

include("representation/cargo.jl")
include("representation/deck.jl")
include("representation/slot.jl")
include("representation/vessel.jl")
include("representation/instance.jl")
include("representation/problem.jl")

include("model/base_model.jl")
include("model/hazardous.jl")
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

#include("loadmaster/types.jl")
#include("loadmaster/api.jl")


# Makes these functions available outside the module
# load data and solve model function - Deterministic
export load_data, create_model, solve_model, extract_solution, plot_solution
# Stochastic functions
export create_model_stochastic, extract_stochastic_solution, plot_stochastic_solution
export create_stochastic_problem, create_stochastic_model, Monto_Carlo_sampling
export compare_solutions_print, plot_ballast_weight_diff, plot_ballast_cargo_weight_diff
# Struct for stowage problem and solution
export StowageProblem, Solution
# Stochastic structs
export StochasticStowageProblem, SolutionStochastic
# plots
export plot_cargo_weights, plot_cargo_OG

# Functions to connect to LoadMaster
#export LoadMaster, getcurrentvessel, resetvessel, upload_solution, getresults, parse_result, uri, headers, struct_to_json

end # module StowagePlannerStochastic
