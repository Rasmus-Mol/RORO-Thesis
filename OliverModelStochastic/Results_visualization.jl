# Script to plot results
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

# load data from solutions
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