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
# HPC_folder
HPC_folder = "Finlandia_01_04_09_41_38" 
# Load data from HPC
repetitions, scenarios, n_unknown, time_limit = get_HPC_data(HPC_folder)
#repetitions = 1 # didn't finish running
sc = length(scenarios)
n = length(n_unknown)
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous"
Deterministic_problem = load_data(problemname1,problemname2,problemname3)
Deterministic_Solution = get_solution_deterministic("Finlandia_deterministic",
"Deterministic_Solution",HPC_folder)
# soluton arrays
EVP_gen = Array{Any}(undef, repetitions, sc,n)
EVP_gen_fitted = Array{Any}(undef, repetitions, sc,n)
Stochastic_gen = Array{Any}(undef, repetitions, sc,n)
Stochastic_gen_fitted = Array{Any}(undef, repetitions, sc,n)
EVP_boot = Array{Any}(undef, repetitions, sc,n)
EVP_boot_fitted = Array{Any}(undef, repetitions, sc,n)
Stochastic_boot = Array{Any}(undef, repetitions, sc,n)
Stochastic_boot_fitted = Array{Any}(undef, repetitions, sc,n)


for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # Solutions
            # EVP
            foldername = "EVP_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Solution"
            EVP_gen[i,j,k] = get_solution_deterministic(foldername,
            filename,HPC_folder)
            EVP_gen_fitted[i,j,k] = get_solution_deterministic(foldername,
            "Fitted_Solution",HPC_folder)
            foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Solution"
            EVP_boot[i,j,k] = get_solution_deterministic(foldername,
            filename,HPC_folder)
            EVP_boot_fitted[i,j,k] = get_solution_deterministic(foldername,
            "Fitted_Solution",HPC_folder)
            # Stochastic
            foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Solution"
            Stochastic_gen[i,j,k] = get_solution_stochastic(foldername,
            filename,HPC_folder)
            Stochastic_gen_fitted[i,j,k] = get_solution_deterministic(foldername,
            "Fitted_Solution",HPC_folder)
            foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Solution"
            Stochastic_boot[i,j,k] = get_solution_stochastic(foldername,
            filename,HPC_folder)
            Stochastic_boot_fitted[i,j,k] = get_solution_deterministic(foldername,
            "Fitted_Solution",HPC_folder)
        end
    end
end

println("##########################")
println("Total number of models :", repetitions*sc*n)
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            println("############################")
            println("EVP model gen, Cargo loaded: ", EVP_gen[i,j,k].n_cargo_loaded)
            println("EVP model gen after realization, Cargo loaded: ", EVP_gen_fitted[i,j,k].n_cargo_loaded)
            println("Stochastic model gen, Cargo loaded: ", Stochastic_gen[i,j,k].n_cargo_loaded)
            println("Stochastic model gen after realization, Cargo loaded: ", Stochastic_gen_fitted[i,j,k].n_cargo_loaded)
            println("EVP model boot, Cargo loaded: ", EVP_boot[i,j,k].n_cargo_loaded)
            println("EVP model boot after realization, Cargo loaded: ", EVP_boot_fitted[i,j,k].n_cargo_loaded)
            println("Stochastic model boot, Cargo loaded: ", Stochastic_boot[i,j,k].n_cargo_loaded)
            println("Stochastic model boot after realization, Cargo loaded: ", Stochastic_boot_fitted[i,j,k].n_cargo_loaded)
        end
    end
end

# Plot cargo loaded
n_models = 4
r = 1 # choose repetition number
cargo_loaded_EVP_gen = zeros(Int64,sc,n)
cargo_loaded_EVP_gen_fitted = zeros(Int64,sc,n)
cargo_loaded_Stochastic_gen = zeros(Int64,sc,n)
cargo_loaded_Stochastic_gen_fitted = zeros(Int64,sc,n)
cargo_loaded_EVP_boot = zeros(Int64,sc,n)
cargo_loaded_EVP_boot_fitted = zeros(Int64,sc,n)
cargo_loaded_Stochastic_boot = zeros(Int64,sc,n)
cargo_loaded_Stochastic_boot_fitted = zeros(Int64,sc,n)
gap_EVP_gen = zeros(Float64,sc,n)
gap_Stochastic_gen = zeros(Float64,sc,n)
gap_EVP_boot = zeros(Float64,sc,n)
gap_Stochastic_boot = zeros(Float64,sc,n)
for i in 1:sc
    for j in 1:n
        cargo_loaded_EVP_gen[i,j] = Int(EVP_gen[r,i,j].n_cargo_loaded)
        cargo_loaded_EVP_gen_fitted[i,j] = Int(EVP_gen_fitted[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_gen[i,j] = Int(Stochastic_gen[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_gen_fitted[i,j] = Int(Stochastic_gen_fitted[r,i,j].n_cargo_loaded)
        cargo_loaded_EVP_boot[i,j] = Int(EVP_boot[r,i,j].n_cargo_loaded)
        cargo_loaded_EVP_boot_fitted[i,j] = Int(EVP_boot_fitted[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_boot[i,j] = Int(Stochastic_boot[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_boot_fitted[i,j] = Int(Stochastic_boot_fitted[r,i,j].n_cargo_loaded)
        gap_EVP_gen[i,j] = EVP_gen[r,i,j].gap
        gap_Stochastic_gen[i,j] = Stochastic_gen[r,i,j].gap
        gap_EVP_boot[i,j] = EVP_boot[r,i,j].gap
        gap_Stochastic_boot[i,j] = Stochastic_boot[r,i,j].gap
    end
end
# Display loaded cargo
# EVP gen
linew = 1
p = plot(scenarios,cargo_loaded_EVP_gen[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Cargo loaded for different models,\n before realization of cargo weight", label = "EVP-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_EVP_boot[:,end],label = "EVP-Boot", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_gen[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_boot[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.n_cargo_loaded,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig("Plots/Results/CargoLoaded_BeforeRealization.png")
p = plot(scenarios,cargo_loaded_EVP_gen_fitted[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Cargo loaded for different models,\n after realization of cargo weight", label = "EVP-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_EVP_boot_fitted[:,end],label = "EVP-Boot", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_gen_fitted[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_boot_fitted[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.n_cargo_loaded,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig("Plots/Results/CargoLoaded_AfterRealization.png")

# Display Gaps
p = plot(scenarios,gap_EVP_gen[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Gap for different models,\n before realization of cargo weight", label = "EVP-Gen", linewidth =linew)
plot!(scenarios,gap_EVP_boot[:,end],label = "EVP-Boot", linewidth =linew)
plot!(scenarios,gap_Stochastic_gen[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,gap_Stochastic_boot[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.gap,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig("Plots/Results/Gap_BeforeRealization.png")
# Zoomed in
p = plot(scenarios,gap_EVP_gen[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Gap for different models,\n before realization of cargo weight", label = "EVP-Gen", linewidth =linew)
plot!(scenarios,gap_EVP_boot[:,end],label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,gap_Stochastic_gen[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,gap_Stochastic_boot[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.gap,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig("Plots/Results/Gap_BeforeRealization_zoomed.png")

# Number of times problem was unfeasible
count_infeasible = 0
count_infeasible_fitted = 0
count_timelimit_fitted = 0
count_timelimit = 0
problems_infeasible = []
problems_timelimit = []
for i in 1:sc
    for j in 1:n
        # print if not optimal
        if EVP_gen[r,i,j].status != "OPTIMAL"
            if EVP_gen[r,i,j].status == "TIME_LIMIT"
                count_timelimit += 1
                push!(problems_timelimit,EVP_gen[r,i,j])
            else
                count_infeasible += 1
                push!(problems_infeasible,EVP_gen[r,i,j])
            end
            println("##########################")
            println("EVP gen model unfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", EVP_gen[r,i,j].status)
        end
        if EVP_gen_fitted[r,i,j].status != "OPTIMAL"
            if EVP_gen_fitted[r,i,j].status == "TIME_LIMIT"
                count_timelimit_fitted += 1
                push!(problems_timelimit,EVP_gen_fitted[r,i,j])
            else
                count_infeasible_fitted += 1
                push!(problems_infeasible,EVP_gen_fitted[r,i,j])
            end
            println("##########################")
            println("EVP gen fitted model unfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", EVP_gen_fitted[r,i,j].status)
        end
        if Stochastic_gen[r,i,j].status != "OPTIMAL"
            if Stochastic_gen[r,i,j].status == "TIME_LIMIT"
                count_timelimit += 1
                push!(problems_timelimit,Stochastic_gen[r,i,j])
            else
                count_infeasible += 1
                push!(problems_infeasible,Stochastic_gen[r,i,j])
            end
            println("##########################")
            println("Stochastic gen model unfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", Stochastic_gen[r,i,j].status)
        end
        if Stochastic_gen_fitted[r,i,j].status != "OPTIMAL"
            if Stochastic_gen_fitted[r,i,j].status == "TIME_LIMIT"
                count_timelimit_fitted += 1
                push!(problems_timelimit,Stochastic_gen_fitted[r,i,j])
            else
                count_infeasible_fitted += 1
                push!(problems_infeasible,Stochastic_gen_fitted[r,i,j])
            end
            println("##########################")
            println("Stochastic gen fitted model unfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", Stochastic_gen_fitted[r,i,j].status)
        end
        if EVP_boot[r,i,j].status != "OPTIMAL"
            if EVP_boot[r,i,j].status == "TIME_LIMIT"
                count_timelimit += 1
                push!(problems_timelimit,EVP_boot[r,i,j])
            else
                count_infeasible += 1
                push!(problems_infeasible,EVP_boot[r,i,j])
            end
            println("##########################")
            println("EVP boot model unfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", EVP_boot[r,i,j].status)
        end
        if EVP_boot_fitted[r,i,j].status != "OPTIMAL"
            if EVP_boot_fitted[r,i,j].status == "TIME_LIMIT"
                count_timelimit_fitted += 1
                push!(problems_timelimit,EVP_boot_fitted[r,i,j])
            else
                count_infeasible_fitted += 1
                push!(problems_infeasible,EVP_boot_fitted[r,i,j])
            end
            println("##########################")
            println("EVP boot fitted model unfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", EVP_boot_fitted[r,i,j].status)
        end
        if Stochastic_boot[r,i,j].status != "OPTIMAL"
            if Stochastic_boot[r,i,j].status == "TIME_LIMIT"
                count_timelimit += 1
                push!(problems_timelimit,Stochastic_boot[r,i,j])
            else
                count_infeasible += 1
                push!(problems_infeasible,Stochastic_boot[r,i,j])
            end
            println("##########################")
            println("Stochastic boot model unfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", Stochastic_boot[r,i,j].status)
        end
        if Stochastic_boot_fitted[r,i,j].status != "OPTIMAL"
            if Stochastic_boot_fitted[r,i,j].status == "TIME_LIMIT"
                count_timelimit_fitted += 1
                push!(problems_timelimit,Stochastic_boot_fitted[r,i,j])
            else
                count_infeasible_fitted += 1
                push!(problems_infeasible,Stochastic_boot_fitted[r,i,j])
            end
            println("##########################")
            println("Stochastic boot fitted model unfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", Stochastic_boot_fitted[r,i,j].status)
        end
    end
end
println("##########################")
println("Number of times problem was infeasible: ", count_infeasible)
println("Number of times problem was infeasible due to time limit: ", count_timelimit)
println("Number of times problem after realization was infeasible: ", count_infeasible_fitted)
println("Number of times problem after realization was infeasible due to time limit: ", count_timelimit_fitted)

# print problems that have issues
sto_timelimit = []
sto_names = []
EVP_timelimit = []
EVP_names = []
println("##########################")
println("Deterministic total cargo weight: ", Deterministic_Solution.cargo_weight)
for i in 1:length(problems_timelimit)
    println("##########################")
    if typeof(problems_timelimit[i]) == SolutionStochastic
        push!(sto_timelimit, problems_timelimit[i])
        println("Stochastic model")
        println("Average total weight of cargo: ", sum(c.cargo_weight for c in problems_timelimit[i].forces)/length(problems_timelimit[i].forces))
        println("Number of scenarios: ", length(problems_timelimit[i].forces))
    else
        push!(EVP_timelimit, problems_timelimit[i])
        println("EVP model")
        println("Total weight of cargo: ", problems_timelimit[i].cargo_weight)
    end
end
println("##########################")
println("Number of times a EVP-problem reached time: ",length(EVP_timelimit))
println("Number of times a Stochastic-problem reached time: ",length(sto_timelimit))

