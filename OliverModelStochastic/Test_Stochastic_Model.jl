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
# Which instance of HPC results to load - has to be changed manually
HPC_folder = "Finlandia_09_04_09"
# load data
repetitions, scenarios, n_cargo_unknownweight = get_HPC_data(HPC_folder)
# Load problems and solutions
problemname1, problemname2, problemname3 = "finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous"
det_problem = load_data(problemname1,problemname2,problemname3)
det_solution = get_solution_deterministic("Finlandia_deterministic","Deterministic_Solution",HPC_folder)
# HPC data
repetitions, scenarios, n_unknown, time_limit = get_HPC_data(HPC_folder)
sc = length(scenarios)
n = length(n_unknown)
# Problems and models
problems_sto = Array{Any}(undef, repetitions, sc,n)
problems_EVP = Array{Any}(undef, repetitions, sc,n)
models_sto = Array{Any}(undef, repetitions, sc,n)
models_EVP = Array{Any}(undef, repetitions, sc,n)
# Solution
solutions_sto = Array{Any}(undef, repetitions, sc,n)
solutions_EVP = Array{Any}(undef, repetitions, sc,n)
fitted_sol = Array{Any}(undef, repetitions, sc,n)
fitted_sol_EVP = Array{Any}(undef, repetitions, sc,n)
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # EVP
            foldername = "EVP_rep$(i)_sc$(j)_unknown$(k)"
            filename

            # Stochastic

        end
    end
end


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

include("src/model/random_stowage_plan.jl")

@unpack vessel, slots, cargo = det_problem

# choose n first elements and shuffle them
n = 117
cargoc = CargoCollection(shuffle!([cargo.items[i] for i in 1:n]))
splan, not_s = random_stowage_plan(cargoc, slots)
println("Number of cargoes loaded: ", n-length(not_s),"/",n)
#splan2, not_s2 = random_stowage_plan(sort_cargocollection(cargoc), slots)
# New random cargocollection
ntypes = [length(filter(x -> x.cargo_type_id == i, cargoc.items)) for i in 1:4]
rcargoc = random_cargocollection(ntypes)
splan_r, not_s_r = random_stowage_plan(rcargoc, slots)
println("Number of cargoes loaded: ", n-length(not_s_r),"/",n)
#splan_r2, not_s_r2 = random_stowage_plan(sort_cargocollection(rcargoc), slots)

# Original solution
plot_solution(det_solution) 
# Random stowage plan
plot_solution_random_plan(splan, cargoc,slots)
#plot_solution_random_plan(splan2, sort_cargocollection(cargoc),slots)
# Another random stowage plan
plot_solution_random_plan(splan_r, rcargoc,slots)
#plot_solution_random_plan(splan_r2, sort_cargocollection(rcargoc),slots)  

# Test how many cargoes need to be moved to make feasible plan
mo = create_random_stowageplan_model(splan, not_s, cargoc, vessel, slots, false)
set_silent(mo) # removes terminal output
set_time_limit_sec(mo, 60)
optimize!(mo)
pro = StowageProblem(vessel = vessel,slots = slots,cargo = cargoc,name = det_problem.name)
sol = extract_solution(pro,mo)
sol_plan = sol.cs
y = value.(mo[:y])
println("Moves: ", sum(y))
not_loaded = [i[1] for i in findall(x -> x<0.5, sum(sol_plan, dims = 2))]
not_loaded_random_plan = [i[1] for i in findall(x -> x<0.5, sum(splan, dims = 2))]
println("Not loaded after optimization: ", length(not_loaded))
println("Not loaded before optimization: ", length(not_s))
println("Difference: ", setdiff(not_loaded, not_loaded_random_plan))

# plot_solution(sol) Can't use this function to show solution
plot_solution_random_plan(sol_plan, cargoc,slots) 
