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

# Creates Deterministic problem and model
problem_det = load_data("finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous")
model_det = create_model(problem_det)
set_silent(model_det) # removes terminal output
set_time_limit_sec(model_det, 60 * 5)
optimize!(model_det)
solution_det = extract_solution(problem_det, model_det)
println("Deterministic model time: ", solution_det.time)
println("Deterministic model cargo loaded: ", solution_det.n_cargo_loaded)
println("Deterministic model status: ", solution_det.status)

# Creates Stochastic problem and model
# parameters
scenarios = [10,20,50]
n_cargo_unknownweight = [10,50,length(problem_det.cargo)] # all cargo weights are unknown
sc = scenarios[2]
n_c_un = n_cargo_unknownweight[2]

time_limit = 60 * 5
problems_gen = []
problems_monte = []
models_gen = []
models_monte = []
sol_gen = []
sol_monte = []
cs_gen = []
cs_monte = []
fitted_sol_gen = []
fitted_sol_monte = []

repetitions = 5
for i in 1:repetitions
    println("Iteration: ", i)
    pro_sto = create_stochastic_problem(problem_det, sc, n_c_un,[])
    push!(problems_gen, pro_sto)
    model_sto = create_model_stochastic(pro_sto)
    push!(models_gen, model_sto)
    set_silent(model_sto) # removes terminal output
    set_time_limit_sec(model_sto, time_limit) # 5 minutes to solve model
    optimize!(model_sto)
    sol_sto = extract_stochastic_solution(pro_sto, model_sto)
    push!(sol_gen, sol_sto)
    cs_sto = sol_sto.cs
    push!(cs_gen, cs_sto)
    second_stage_m = second_stage_model(cs_sto, problem_det)
    set_silent(second_stage_m) # removes terminal output
    set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
    optimize!(second_stage_m)
    fitted_sol = get_solution_second_stage(problem_det, second_stage_m, sol_sto)
    push!(fitted_sol_gen, fitted_sol)

    pro_sto = create_stochastic_problem(problem_det, sc, n_c_un,[],Monto_Carlo_sampling)
    push!(problems_monte, pro_sto)
    model_sto = create_model_stochastic(pro_sto)
    push!(models_monte, model_sto)
    set_silent(model_sto) # removes terminal output
    set_time_limit_sec(model_sto, time_limit) # 5 minutes to solve model
    optimize!(model_sto)
    sol_sto = extract_stochastic_solution(pro_sto, model_sto)
    push!(sol_monte, sol_sto)
    cs_sto = sol_sto.cs
    push!(cs_monte, cs_sto)
    second_stage_m = second_stage_model(cs_sto, problem_det)
    set_silent(second_stage_m) # removes terminal output
    set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
    optimize!(second_stage_m)
    fitted_sol = get_solution_second_stage(problem_det, second_stage_m, sol_sto)
    push!(fitted_sol_monte, fitted_sol)
end


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

# Total cargo weight
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
        label="Repetition $i", xlabel="Scenarios", ylabel="Total argo weight difference (t)", 
    title="Total cargo weight difference,\n scenarios: $sc, unknown weight: $n_c_un")
    else
        plot!(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_gen[i].cargo.items],label="Repetition $i")
    end
end
savefig("Results_Plots/totalweightdiffAll_gen_sc_$(sc)_n_c_un_$(n_c_un).png")
for i in 1:repetitions
    if i == 1
        p = plot(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_monte[i].cargo.items],
        label="Repetition $i", xlabel="Scenarios", ylabel="Total argo weight difference (t)", 
    title="Total cargo weight difference,\n scenarios: $sc, unknown weight: $n_c_un")
    else
        plot!(1:sc, [item.total_weight-problem_det.cargo.total_weight for item in problems_monte[i].cargo.items],label="Repetition $i")
    end
end
savefig("Results_Plots/totalweightdiffAll_monte_sc_$(sc)_n_c_un_$(n_c_un).png")



