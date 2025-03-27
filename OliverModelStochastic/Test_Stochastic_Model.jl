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


include("src/representation/cargo.jl")
include("src/representation/deck.jl")
include("src/representation/slot.jl")
include("src/representation/vessel.jl")
include("src/representation/instance.jl")
include("src/representation/problem.jl")
# New: Stochastic 
include("src/representation/CargoScenarios.jl")
include("src/plots/weight_plots.jl")
include("src/model/base_model.jl")
include("src/model/stability.jl")
include("src/model/base_stochastic_model.jl")
include("src/model/stability_stochastic.jl")
include("src/model/second_stage_model.jl")
include("src/utils/helpers.jl")
include("src/solution.jl")
include("src/CompareSolutions.jl")
include("src/plots/solution.jl")
include("src/utils/SaveData.jl")

################################################################
# Load problems and solutions
problemname1, problemname2, problemname3 = "finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous"
det_problem = load_data(problemname1,problemname2,problemname3)
det_solution = get_solution_deterministic("Finlandia_deterministic","Deterministic_Solution")
sto_problem = []
sto_solution = []
EVP_problem = []
EVP_solution = []

get_deterministic_problem("",filename::String,problemname1,problemname2,problemname3)
get_stochastic_problem(foldername::String,filename::String,problemname1,problemname2,problemname3)
get_solution_deterministic(foldername::String,filename::String)
################################################################

# plot
plot(1:repetitions, [sol.time for sol in sol_gen], label="Generalized scenario", xlabel="Repetitions", ylabel="Time (s)", title="Time to solve model,\n scenarios: $sc, unknown weight: $n_c_un")
savefig("Results_Plots/time_gen_sc_$(sc)_n_c_un_$(n_c_un).png")
plot(1:repetitions, [sol.time for sol in sol_monte], label="Monte Carlo scenario", xlabel="Repetitions", ylabel="Time (s)", title="Time to solve model,\n scenarios: $sc, unknown weight: $n_c_un")
savefig("Results_Plots/time_monte_sc_$(sc)_n_c_un_$(n_c_un).png")
plot(1:repetitions, [sol.n_cargo_loaded for sol in fitted_sol_gen], label="Generalized scenario", xlabel="Repetitions", ylabel="Cargo loaded", title="Cargo loaded,\n scenarios: $sc, unknown weight: $n_c_un")
savefig("Results_Plots/cargoloaded_gen_sc_$(sc)_n_c_un_$(n_c_un).png")
plot(1:repetitions, [sol.n_cargo_loaded for sol in fitted_sol_monte], label="Monte Carlo scenario", xlabel="Repetitions", ylabel="Cargo loaded", title="Cargo loaded,\n scenarios: $sc, unknown weight: $n_c_un")
savefig("Results_Plots/cargoloaded_monte_sc_$(sc)_n_c_un_$(n_c_un).png")
plot(1:repetitions, [fitted_sol.ballast_weight-solution_det.ballast_weight for fitted_sol in fitted_sol_gen], label="Generalized scenario", xlabel="Repetitions", ylabel="Ballast weight (t)", title="Difference in ballast weight,\n scenarios: $sc, unknown weight: $n_c_un")
savefig("Results_Plots/ballastdif_gen_sc_$(sc)_n_c_un_$(n_c_un).png")
plot(1:repetitions, [fitted_sol.ballast_weight-solution_det.ballast_weight for fitted_sol in fitted_sol_monte], label="Monte Carlo scenario", xlabel="Repetitions", ylabel="Ballast weight (t)", title="Difference in ballast weight,\n scenarios: $sc, unknown weight: $n_c_un")
savefig("Results_Plots/ballastdif_monte_sc_$(sc)_n_c_un_$(n_c_un).png")

# Total cargo weight compared to deterministic
for i in 1:repetitions
    plot(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_gen[i].cargo.items], 
    label="Repetition $i", xlabel="Scenarios", ylabel="Total cargo weight difference (t)", 
    title="Total cargo weight difference,\n scenarios: $sc, unknown weight: $n_c_un, rep: $i")
    for j in 1:sc
        temp = problems_gen[i].cargo.items[j].total_weight-problem_det.cargo.total_weight
        annotate!(j, temp+20 , text(round(temp,digits = 1),7))
    end
    savefig("Results_Plots/totalweightdiff_gen_sc_$(sc)_n_c_un_$(n_c_un)_rep_$(i).png")

    plot(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_monte[i].cargo.items], 
    label="Repetition $i", xlabel="Scenarios", ylabel="Total cargo weight difference (t)", 
    title="Total cargo weight difference,\n scenarios: $sc, unknown weight: $n_c_un, rep: $i")
    for j in 1:sc
        temp = problems_monte[i].cargo.items[j].total_weight-problem_det.cargo.total_weight
        annotate!(j, temp+10 , text(round(temp,digits = 1),7))
    end
    savefig("Results_Plots/totalweightdiff_monte_sc_$(sc)_n_c_un_$(n_c_un)_rep_$(i).png")
end
# all repetitions in same plot
for i in 1:repetitions
    if i == 1
        p = plot(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_gen[i].cargo.items],
        label="Repetition $i", xlabel="Scenarios", ylabel="Total cargo weight difference (t)", 
    title="Total cargo weight difference,\n scenarios: $sc, unknown weight: $n_c_un")
    else
        plot!(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_gen[i].cargo.items],label="Repetition $i")
    end
end
savefig("Results_Plots/totalweightdiffAll_gen_sc_$(sc)_n_c_un_$(n_c_un).png")
for i in 1:repetitions
    if i == 1
        p = plot(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_monte[i].cargo.items],
        label="Repetition $i", xlabel="Scenarios", ylabel="Total cargo weight difference (t)", 
    title="Total cargo weight difference,\n scenarios: $sc, unknown weight: $n_c_un")
    else
        plot!(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_monte[i].cargo.items],label="Repetition $i")
    end
end
savefig("Results_Plots/totalweightdiffAll_monte_sc_$(sc)_n_c_un_$(n_c_un).png")

################################################################

temp = problems_gen[1]
EVP = expected_value_problem(temp)
