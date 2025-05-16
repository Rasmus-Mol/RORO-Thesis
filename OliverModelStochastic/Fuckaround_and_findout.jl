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
include("src/representation/ScenarioReduction.jl")

# Do stuff 
test_instance_hol = Hollandia_test[5]
test_instance_fin = Finlandia_test[6]

Deterministic_problem_hol = load_data("hollandia",test_instance_hol,"hazardous")
#Deterministic_problem_fin = load_data("finlandia","no_cars_medium_100_haz_eq_0.1","hazardous")
Deterministic_problem_fin = load_data("finlandia",test_instance_fin,"hazardous")

pro_hol = create_stochastic_problem(Deterministic_problem_hol, 40, length(Deterministic_problem_hol.cargo), []) 
pro_fin = create_stochastic_problem(Deterministic_problem_fin, 10, length(Deterministic_problem_fin.cargo), []) 


vessel = pro_fin.vessel
vessel.ballast_tanks[1].max_vol
test = calculate_vcg_slopes(vessel)
cargo = Deterministic_problem_fin.cargo

# Hollandia
model_hol = create_model(Deterministic_problem_hol)
set_silent(model_hol) # removes terminal output
set_time_limit_sec(model_hol, 60 * 15) # 5 minutes to solve model
optimize!(model_hol)
model_stochastic_hol = create_model_stochastic(pro_hol)
set_silent(model_stochastic_hol) # removes terminal output
set_time_limit_sec(model_stochastic_hol, 60 * 2) # 5 minutes to solve model
#set_time_limit_sec(model_stochastic, 60 * 15) # 10 minutes to solve model
optimize!(model_stochastic_hol)
sol_hol = extract_stochastic_solution(pro_hol,model_stochastic_hol)
folder1 = "Temp1_delete"
folder2 = "Temp2_delete"
write_solution_stochastic(sol_hol,folder2,"Stochastic_Solution",folder1)
cs_sol = sol_hol.cs
second_stage_m = second_stage_model(cs_sol, Deterministic_problem_hol)
set_silent(second_stage_m) # removes terminal output
set_time_limit_sec(second_stage_m, 60) 
optimize!(second_stage_m)
fitted_sol = get_solution_second_stage_stochastic(Deterministic_problem_hol, second_stage_m, sol_hol)
write_solution(fitted_sol,folder1,"Fitted_Solution",folder2)

# Finlandia 
model_fin = create_model(Deterministic_problem_fin)
set_silent(model_fin) # removes terminal output
set_time_limit_sec(model_fin, 60 * 5) # 5 minutes to solve model
optimize!(model_fin)
solution_det = extract_solution(Deterministic_problem_fin, model_fin)

# med ingen ændringer
# 89 sek med korrekt vcg slopes
# 41 sek med korrekt z_min og z_max
# med begge ændringer

model_stochastic_fin = create_model_stochastic(pro_fin)
set_silent(model_stochastic_fin) # removes terminal output
set_time_limit_sec(model_stochastic_fin, 60 * 5) # 5 minutes to solve model
#set_time_limit_sec(model_stochastic, 60 * 15) # 10 minutes to solve model
optimize!(model_stochastic_fin)
sol = extract_stochastic_solution(pro_fin,model_stochastic_fin)
cs_sol = sol.cs

# Solve second stage when we know unknown weights
second_stage_m = second_stage_model(cs_sol, problem_det)
set_silent(second_stage_m) # removes terminal output
set_time_limit_sec(second_stage_m, time_limit) 
optimize!(second_stage_m)

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

for i in 1:length(Hollandia_test)
    test_instance_hol = Hollandia_test[i]
    test_instance_fin = Finlandia_test[i]
Deterministic_problem_hol = load_data("hollandia",test_instance_hol,"hazardous")
Deterministic_problem_fin = load_data("finlandia",test_instance_fin,"hazardous")

    println("Finlandia test $(i)")
    println("Total weight: ", Deterministic_problem_fin.cargo.total_weight)
    println("number of cars:", length(filter(x -> x.cargo_type_id == 2, Deterministic_problem_fin.cargo)))
    println("number of trucks:", length(filter(x -> x.cargo_type_id == 1, Deterministic_problem_fin.cargo)))
    println("number of Secu-boxes:", length(filter(x -> x.cargo_type_id == 4, Deterministic_problem_fin.cargo)))
    println("number of heavy machinery:", length(filter(x -> x.cargo_type_id == 3, Deterministic_problem_fin.cargo)))
    
    println("Hollandia test $(i)")
    println("Total weight: ", Deterministic_problem_hol.cargo.total_weight)
    println("number of cars:", length(filter(x -> x.cargo_type_id == 2, Deterministic_problem_hol.cargo)))
    println("number of trucks:", length(filter(x -> x.cargo_type_id == 1, Deterministic_problem_hol.cargo)))
    println("number of Secu-boxes:", length(filter(x -> x.cargo_type_id == 4, Deterministic_problem_hol.cargo)))
    println("number of heavy machinery:", length(filter(x -> x.cargo_type_id == 3, Deterministic_problem_hol.cargo))) 
end