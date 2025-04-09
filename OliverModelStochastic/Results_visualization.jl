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

# model
include("src/representation/CargoScenarios.jl")
include("src/plots/weight_plots.jl")
include("src/model/base_model.jl")
include("src/model/stability.jl")
include("src/model/base_stochastic_model.jl")
include("src/model/stability_stochastic.jl")
include("src/model/second_stage_model.jl")
include("src/utils/helpers.jl")

# load data from solutions
# HPC_folder
#HPC_folder = "Finlandia_03_04_13"
#plot_folder = "Plots/Results/day_03_04_version_bootstrap_1/"
#HPC_folder = "Finlandia_07_04_14"
#plot_folder = "Plots/Results/day_07_04_version_bootstrap_1_first_results/"
HPC_folder = "Finlandia_07_04_15"
plot_folder = "Plots/Results/day_07_04_version_bootstrap_1_second_results/"


# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end
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

test = Stochastic_boot_fitted[1,3,1]
test.status == "INFEASIBLE"

println("##########################")
println("Total number of models :", repetitions*sc*n)
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            println("############################")
            println("Number of scnearios: ", scenarios[j])
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
savefig(plot_folder*"CargoLoaded_BeforeRealization.png")
p = plot(scenarios,cargo_loaded_EVP_gen_fitted[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Cargo loaded for different models,\n after realization of cargo weight", label = "EVP-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_EVP_boot_fitted[:,end],label = "EVP-Boot", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_gen_fitted[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_boot_fitted[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.n_cargo_loaded,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"CargoLoaded_AfterRealization.png")

# Display Gaps
p = plot(scenarios,gap_EVP_gen[:,end] .*100,xlabel="Scenarios",ylabel="Gap (%)",
title="Gap for different models,\n before realization of cargo weight", label = "EVP-Gen", linewidth =linew)
plot!(scenarios,gap_EVP_boot[:,end] .*100,label = "EVP-Boot", linewidth =linew)
plot!(scenarios,gap_Stochastic_gen[:,end] .* 100,label = "Sto-Gen", linewidth =linew)
plot!(scenarios,gap_Stochastic_boot[:,end] .* 100,label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.gap .*100,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"Gap_BeforeRealization.png")
# Zoomed in
p = plot(scenarios,gap_EVP_gen[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Gap for different models,\n before realization of cargo weight", label = "EVP-Gen", linewidth =linew)
plot!(scenarios,gap_EVP_boot[:,end],label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,gap_Stochastic_gen[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,gap_Stochastic_boot[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.gap,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"Gap_BeforeRealization_zoomed.png")

# Number of times problem was unfeasible
EVP_gen_inf = []
EVP_gen_fitted_inf = []
Stochastic_gen_inf = []
Stochastic_gen_fitted_inf = []
EVP_boot_inf = []
EVP_boot_fitted_inf = []
Stochastic_boot_inf = []
Stochastic_boot_fitted_inf = []

for i in 1:sc
    for j in 1:n
        # print if not optimal
        if EVP_gen[r,i,j].status != "OPTIMAL"
            push!(EVP_gen_inf,EVP_gen[r,i,j])
            println("##########################")
            println("EVP gen model infeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", EVP_gen[r,i,j].status)
        end
        if EVP_gen_fitted[r,i,j].status != "OPTIMAL"
            push!(EVP_gen_fitted_inf,EVP_gen_fitted[r,i,j])
            println("##########################")
            println("EVP gen fitted model jnfeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", EVP_gen_fitted[r,i,j].status)
        end
        if Stochastic_gen[r,i,j].status != "OPTIMAL"
            push!(Stochastic_gen_inf,Stochastic_gen[r,i,j])
            println("##########################")
            println("Stochastic gen model infeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", Stochastic_gen[r,i,j].status)
        end
        if Stochastic_gen_fitted[r,i,j].status != "OPTIMAL"
            push!(Stochastic_gen_fitted_inf,Stochastic_gen_fitted[r,i,j])
            println("##########################")
            println("Stochastic gen fitted model infeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", Stochastic_gen_fitted[r,i,j].status)
        end
        if EVP_boot[r,i,j].status != "OPTIMAL"
            push!(EVP_boot_inf,EVP_boot[r,i,j])
            println("##########################")
            println("EVP boot model infeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", EVP_boot[r,i,j].status)
        end
        if EVP_boot_fitted[r,i,j].status != "OPTIMAL"
            push!(EVP_boot_fitted_inf,EVP_boot_fitted[r,i,j])
            println("##########################")
            println("EVP boot fitted model infeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", EVP_boot_fitted[r,i,j].status)
        end
        if Stochastic_boot[r,i,j].status != "OPTIMAL"
            push!(Stochastic_boot_inf,Stochastic_boot[r,i,j])
            println("##########################")
            println("Stochastic boot model infeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", Stochastic_boot[r,i,j].status)
        end
        if Stochastic_boot_fitted[r,i,j].status != "OPTIMAL"
            push!(Stochastic_boot_fitted_inf,Stochastic_boot_fitted[r,i,j])
            println("##########################")
            println("Stochastic boot fitted model infeasible for repetition $(r), scenario $(i), unknown cargo $(j)")
            println("Status: ", Stochastic_boot_fitted[r,i,j].status)
        end
    end
end
println("##########################")
println("Number of different parameters: ", sc*n)
println("Stochastic models + EVP models: ", length(EVP_gen)+length(EVP_boot)+length(Stochastic_gen)+length(Stochastic_boot))
println("##########################")
# EVP models - normally none of them should be infeasible
for i in 1:length(EVP_gen_inf)
    println("EVP-gen Status: ", EVP_gen_inf[i].status)
    println("Index: ", i)
end
println("##########################")
for i in 1:length(EVP_boot_inf)
    println("EVP-boot Status: ", EVP_boot_inf[i].status)
    println("Index: ", i)
end
# Stochastic models 
println("##########################")
println("Stochastic gen model:")
for i in 1:length(Stochastic_gen_inf)
    println("Stochastic gen model infeasible for number of scenarios: $(length(Stochastic_gen_inf[i].forces))")
    println("Status: ", Stochastic_gen_inf[i].status)
    if Stochastic_gen_inf[i].status == "TIME_LIMIT"
        println("Cargo loaded: ", Stochastic_gen_inf[i].n_cargo_loaded)
    else # if infeasible
        for j in 1:lengt(h(Stochastic_gen_inf[i].forces))
            println("Weight in scenarios", Stochastic_gen_inf[i].forces[j].cargo_weight)
        end
    end
end
println("##########################")
println("Stochastic Boot model:")
for i in 1:length(Stochastic_boot_inf)
    println("Stochastic boot model infeasible for number of scenarios: $(length(Stochastic_boot_inf[i].forces))")
    println("Status: ", Stochastic_boot_inf[i].status)
    if Stochastic_boot_inf[i].status == "TIME_LIMIT"
        println("Cargo loaded: ", Stochastic_boot_inf[i].n_cargo_loaded)
    else # if infeasible
        for j in 1:lengt(h(Stochastic_boot_inf[i].forces))
            println("Weight in scenarios", Stochastic_boot_inf[i].forces[j].cargo_weight)
        end
    end
end
# After the realization: The recourse Model
println("##########################")
# EVP models - normally none of them should be infeasible
println("EVP-fitted gen model:")
for i in 1:length(EVP_gen_fitted_inf)
    println("Status: ", EVP_gen_fitted_inf[i].status)
end
println("##########################")
println("EVP-fitted boot model:")
for i in 1:length(EVP_boot_fitted_inf)
    println("Status: ", EVP_boot_fitted_inf[i].status)
    if EVP_boot_fitted_inf[i].status != "TIME_LIMIT"
        # TODO Do something to check why this is the case
    end
end
# Stochastic models 
println("##########################")
println("Stochastic-fitted gen model:")
for i in 1:length(Stochastic_gen_fitted_inf)
    println("Stochastic gen model infeasible.")
    println("Status: ", Stochastic_gen_fitted_inf[i].status)
    if Stochastic_gen_fitted_inf[i].status != "TIME_LIMIT" # infeasible
        # TODO Do something to check why this is the case
    end
end
println("##########################")
println("Stochastic-fitted boot model:")
for i in 1:length(Stochastic_boot_fitted_inf)
    println("Stochastic boot model infeasible.")
    println("Status: ", Stochastic_boot_fitted_inf[i].status)
    if Stochastic_boot_fitted_inf[i].status != "TIME_LIMIT"
        # TODO Do something to check why this is the case
    end
end

##################################3
# Check model with slack variables

# Stochastic Boot
slack_variables_boot = []
for i in 1:sc
    foldername = "Stochastic_Bootstrap1_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        # Get placement from stochastic problem
        temp = get_solution_stochastic(foldername,"Stochastic_Solution",HPC_folder)
        cs = temp.cs
        model = second_stage_model_slack(cs, Deterministic_problem)
        set_silent(model) # removes terminal output
        set_time_limit_sec(model, 60 * 5) # 5 minutes to start with
        optimize!(model)
        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("slack model status: ", status)
        else
            push!(slack_variables_boot,[value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax])
            ])
        end
    end
end
# Slack variables
println("##########################")
println("Slack variables for Stochastic Boot")
for i in 1:length(slack_variables_boot)
    for j in 1:length(slack_variables_boot[i])
        if any(slack_variables_boot[i][j] .!= 0.0)
            if j == 1
                println("In model $(i), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_boot[i][j])
            else
                println("Model $(i) and slack variable index $(j) are not zero")
                println(slack_variables_boot[i][j])
            end
        end
    end
end

# EVP Boot
slack_variables_EVP_boot = []
for i in 1:sc
    foldername = "EVP_Bootstrap1_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        # Get placement from stochastic problem
        temp = get_solution_deterministic(foldername,"EVP_Solution",HPC_folder)
        cs = temp.cs
        model = second_stage_model_slack(cs, Deterministic_problem)
        set_silent(model) # removes terminal output
        set_time_limit_sec(model, 60 * 5) # 5 minutes to start with
        optimize!(model)
        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("slack model status: ", status)
        else
            push!(slack_variables_EVP_boot,[value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax])
            ])
        end
    end
end
# Slack variables
println("##########################")
println("Slack variables for EVP Boot")
for i in 1:length(slack_variables_EVP_boot)
    for j in 1:length(slack_variables_EVP_boot[i])
        if any(slack_variables_EVP_boot[i][j] .!= 0.0)
            if j == 1
                println("In model $(i), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_EVP_boot[i][j])
            else
                println("Model $(i) and slack variable index $(j) are not zero")
                println(slack_variables_EVP_boot[i][j])
            end
        end
    end
end
# Stochastic Gen 
slack_variables_gen = []
for i in 1:sc
    foldername = "Stochastic_random_sampling_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    if (temp.status == "INFEASIBLE") # solve slack model
        # Get placement from stochastic problem
        temp = get_solution_stochastic(foldername,"Stochastic_Solution",HPC_folder)
        cs = temp.cs
        model = second_stage_model_slack(cs, Deterministic_problem)
        set_silent(model) # removes terminal output
        set_time_limit_sec(model, 60 * 5) # 5 minutes to start with
        optimize!(model)
        status = termination_status(model)
        println(typeof(status))
        if status != MOI.OPTIMAL
            println("slack model status: ", status)
        else
            push!(slack_variables_gen,[value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax])
            ])
        end
    end
end
# Slack variables
println("##########################")
for i in 1:length(slack_variables_gen)
    for j in 1:length(slack_variables_gen[i])
        if any(slack_variables_gen[i][j] .!= 0.0)
            if j == 1
                println("In model $(i), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_gen[i][j])
            else
                println("Model $(i) and slack variable index $(j) are not zero")
                println(slack_variables_gen[i][j])
            end
        end
    end
end
# EVP Gen 
slack_variables_EVP_gen = []
for i in 1:sc
    foldername = "EVP_random_sampling_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        # Get placement from stochastic problem
        temp = get_solution_deterministic(foldername,"EVP_Solution",HPC_folder)
        cs = temp.cs
        model = second_stage_model_slack(cs, Deterministic_problem)
        set_silent(model) # removes terminal output
        set_time_limit_sec(model, 60 * 5) # 5 minutes to start with
        optimize!(model)
        status = termination_status(model)
        println(typeof(status))
        if status != MOI.OPTIMAL
            println("slack model status: ", status)
        else
            push!(slack_variables_EVP_gen,[value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax])
            ])
        end
    end
end
# Slack variables
println("##########################")
for i in 1:length(slack_variables_EVP_gen)
    for j in 1:length(slack_variables_EVP_gen[i])
        if any(slack_variables_EVP_gen[i][j] .!= 0.0)
            if j == 1
                println("In model $(i), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_EVP_gen[i][j])
            else
                println("Model $(i) and slack variable index $(j) are not zero")
                println(slack_variables_EVP_gen[i][j])
            end
        end
    end
end
