# Script to fuck around in and find out
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

# model
include("src/representation/CargoScenarios.jl")
include("src/plots/weight_plots.jl")
include("src/model/base_model.jl")
include("src/model/stability.jl")
include("src/model/base_stochastic_model.jl")
include("src/model/stability_stochastic.jl")
include("src/model/second_stage_model.jl")
include("src/utils/helpers.jl")
include("src/utils/test_instances.jl")

# Do stuff 
Deterministic_problem_hol = load_data("hollandia","mixed_light_100","hazardous")
Deterministic_problem_fin = load_data("finlandia","no_cars_medium_100","hazardous")
length(Deterministic_problem_fin.slots)

@unpack vessel, slots, cargo = Deterministic_problem_hol
vessel_long = Vessel("finlandia")
vessel_short = simplify_vessel(vessel_long, target_points=50)

c_length = []
for i in 1:length(Finlandia_test)
    #println("Test instance: ", Finlandia_test[i])
    Deterministic_problem = load_data("finlandia", Finlandia_test[i], "hazardous")
    push!(c_length, length(Deterministic_problem.cargo))
end
println("Cargo lengths: ", c_length)
findall(x -> x == 117, c_length)
Finlandia_test[27]
Finlandia_test[31]
Finlandia_test[35]

lcg_all = [x.lcg for x in slots]
minimum(lcg_all)
maximum(lcg_all)

#=
for i in 1:length(vessel.decks)
    println("Deck $(i) weight limit: ", vessel.decks[i].weight_limit)
end

deck1 = filter(x -> x.deck_id == 1, slots)
slots1 = [filter(x->x.cargo_type_id == i, deck1) for i in 1:4]
println("Deck 1")
println("trucks: ", length(slots1[1]))
println("cars: ", length(slots1[2]))
println("Heavy machinery: ", length(slots1[3]))
println("Secu: ", length(slots1[4]))

deck2 = filter(x -> x.deck_id == 2, slots)
slots2 = [filter(x->x.cargo_type_id == i, deck2) for i in 1:4]
println("Deck 2")
println("trucks: ", length(slots2[1]))
println("cars: ", length(slots2[2]))
println("Heavy machinery: ", length(slots2[3]))
println("Secu: ", length(slots2[4]))

deck3 = filter(x -> x.deck_id == 3, slots)
slots3 = [filter(x->x.cargo_type_id == i, deck3) for i in 1:4]
println("Deck 3")
println("trucks: ", length(slots3[1]))
println("cars: ", length(slots3[2]))
println("Heavy machinery: ", length(slots3[3]))
println("Secu: ", length(slots3[4]))

deck4 = filter(x -> x.deck_id == 4, slots)
slots4 = [filter(x->x.cargo_type_id == i, deck4) for i in 1:4]
println("Deck 4")
println("trucks: ", length(slots4[1]))
println("cars: ", length(slots4[2]))
println("Heavy machinery: ", length(slots4[3]))
println("Secu: ", length(slots4[4]))

deck5 = filter(x -> x.deck_id == 5, slots)
slots5 = [filter(x -> x.cargo_type_id == i, deck5) for i in 1:4]
println("Deck 5")
println("trucks: ", length(slots5[1]))
println("cars: ", length(slots5[2]))
println("Heavy machinery: ", length(slots5[3]))
println("Secu: ", length(slots5[4]))
=#
