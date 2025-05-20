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

include("src/utils/test_instances.jl")

# load data from solutions
# HPC_folder
test_instance = Finlandia_test[4]
HPC_folder = "Finlandia_"*test_instance*"_15_05_17"
plot_folder = "Plots/Results/Finlandia_"*test_instance*"/"
# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end
# Load data from HPC
repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
println("Extra info about test: ", note)
#repetitions = 1 # didn't finish running
sc = length(scenarios)
n = length(n_unknown)
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", test_instance, "hazardous"
Deterministic_problem = load_data(problemname1,problemname2,problemname3)
Deterministic_Solution = get_solution_deterministic("Finlandia_deterministic",
"Deterministic_Solution",HPC_folder)
# soluton arrays
#EVP_gen = Array{Any}(undef, repetitions, sc,n)
#EVP_gen_fitted = Array{Any}(undef, repetitions, sc,n)
Stochastic_gen = Array{Any}(undef, repetitions, sc,n)
Stochastic_gen_fitted = Array{Any}(undef, repetitions, sc,n)
#EVP_boot = Array{Any}(undef, repetitions, sc,n)
#EVP_boot_fitted = Array{Any}(undef, repetitions, sc,n)
Stochastic_boot = Array{Any}(undef, repetitions, sc,n)
Stochastic_boot_fitted = Array{Any}(undef, repetitions, sc,n)


for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # Solutions
            # EVP
            #foldername = "EVP_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            #filename = "EVP_Solution"
            #EVP_gen[i,j,k] = get_solution_deterministic(foldername,
            #filename,HPC_folder)
            #EVP_gen_fitted[i,j,k] = get_solution_deterministic(foldername,
            #"Fitted_Solution",HPC_folder)
            #foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            #filename = "EVP_Solution"
            #EVP_boot[i,j,k] = get_solution_deterministic(foldername,
            #filename,HPC_folder)
            #EVP_boot_fitted[i,j,k] = get_solution_deterministic(foldername,
            #"Fitted_Solution",HPC_folder)
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

println("Cargo in problem: ", length(Deterministic_problem.cargo))
println("Deterministic solution packs: ",Deterministic_Solution.n_cargo_loaded)

println("##########################")
println("Total number of models :", repetitions*sc*n)
println("Deterministic solution, Cargo loaded: ", Deterministic_Solution.n_cargo_loaded)
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            println("############################")
            println("Number of scnearios: ", scenarios[j])
            #println("EVP model gen, Cargo loaded: ", EVP_gen[i,j,k].n_cargo_loaded)
            #println("EVP model gen after realization, Cargo loaded: ", EVP_gen_fitted[i,j,k].n_cargo_loaded)
            println("Stochastic model gen, Cargo loaded: ", Stochastic_gen[i,j,k].n_cargo_loaded)
            println("Stochastic model gen after realization, Cargo loaded: ", Stochastic_gen_fitted[i,j,k].n_cargo_loaded)
            #println("EVP model boot, Cargo loaded: ", EVP_boot[i,j,k].n_cargo_loaded)
            #println("EVP model boot after realization, Cargo loaded: ", EVP_boot_fitted[i,j,k].n_cargo_loaded)
            println("Stochastic model boot, Cargo loaded: ", Stochastic_boot[i,j,k].n_cargo_loaded)
            println("Stochastic model boot after realization, Cargo loaded: ", Stochastic_boot_fitted[i,j,k].n_cargo_loaded)
        end
    end
end
println("##########################")
for i in 1:sc
    println("Average number of cargo loaded, scenarios = $(scenarios[i]) - gen: ",
    sum([Stochastic_gen[j,i,end].n_cargo_loaded for j in 1:repetitions])/repetitions)
    println("Average number of cargo loaded, scenarios = $(scenarios[i]) - boot: ",
    sum([Stochastic_boot[j,i,end].n_cargo_loaded for j in 1:repetitions])/repetitions)
end
println("##########################")

# Plot cargo loaded
n_models = 4
r = 1 # choose repetition number
#cargo_loaded_EVP_gen = zeros(Int64,sc,n)
#cargo_loaded_EVP_gen_fitted = zeros(Int64,sc,n)
cargo_loaded_Stochastic_gen = zeros(Int64,sc,n)
cargo_loaded_Stochastic_gen_fitted = zeros(Int64,sc,n)
#cargo_loaded_EVP_boot = zeros(Int64,sc,n)
#cargo_loaded_EVP_boot_fitted = zeros(Int64,sc,n)
cargo_loaded_Stochastic_boot = zeros(Int64,sc,n)
cargo_loaded_Stochastic_boot_fitted = zeros(Int64,sc,n)
#gap_EVP_gen = zeros(Float64,sc,n)
gap_Stochastic_gen = zeros(Float64,sc,n)
#gap_EVP_boot = zeros(Float64,sc,n)
gap_Stochastic_boot = zeros(Float64,sc,n)
for i in 1:sc
    for j in 1:n
        #cargo_loaded_EVP_gen[i,j] = Int(EVP_gen[r,i,j].n_cargo_loaded)
        #cargo_loaded_EVP_gen_fitted[i,j] = Int(EVP_gen_fitted[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_gen[i,j] = Int(Stochastic_gen[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_gen_fitted[i,j] = Int(Stochastic_gen_fitted[r,i,j].n_cargo_loaded)
        #cargo_loaded_EVP_boot[i,j] = Int(EVP_boot[r,i,j].n_cargo_loaded)
        #cargo_loaded_EVP_boot_fitted[i,j] = Int(EVP_boot_fitted[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_boot[i,j] = Int(Stochastic_boot[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_boot_fitted[i,j] = Int(Stochastic_boot_fitted[r,i,j].n_cargo_loaded)
        #gap_EVP_gen[i,j] = EVP_gen[r,i,j].gap
        gap_Stochastic_gen[i,j] = Stochastic_gen[r,i,j].gap
        #gap_EVP_boot[i,j] = EVP_boot[r,i,j].gap
        gap_Stochastic_boot[i,j] = Stochastic_boot[r,i,j].gap
    end
end
# Display loaded cargo
# EVP gen
linew = 1
p = plot(scenarios,cargo_loaded_Stochastic_gen[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Cargo loaded for different models,\n before realization of cargo weight", label = "Sto-Gen", linewidth =linew)
#plot!(scenarios,cargo_loaded_EVP_boot[:,end],label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,cargo_loaded_Stochastic_gen[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_boot[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.n_cargo_loaded,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"CargoLoaded_BeforeRealization.png")
p = plot(scenarios,cargo_loaded_Stochastic_gen_fitted[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Cargo loaded for different models,\n after realization of cargo weight", label = "Sto-Gen", linewidth =linew)
#plot!(scenarios,cargo_loaded_EVP_boot_fitted[:,end],label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,cargo_loaded_Stochastic_gen_fitted[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_boot_fitted[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.n_cargo_loaded,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"CargoLoaded_AfterRealization.png")

# Display Gaps
p = plot(scenarios,gap_Stochastic_gen[:,end] .*100,xlabel="Scenarios",ylabel="Gap (%)",
title="Gap for different models,\n before realization of cargo weight", label = "Sto-Gen", linewidth =linew)
#plot!(scenarios,gap_EVP_boot[:,end] .*100,label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,gap_Stochastic_gen[:,end] .* 100,label = "Sto-Gen", linewidth =linew)
plot!(scenarios,gap_Stochastic_boot[:,end] .* 100,label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.gap .*100,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"Gap_BeforeRealization.png")
# Zoomed in
p = plot(scenarios,gap_Stochastic_boot[:,end].*100,xlabel="Scenarios",ylabel="Gap (%)",
title="Gap for different models,\n before realization of cargo weight", label = "Sto-Boot", linewidth =linew)
#plot!(scenarios,gap_EVP_boot[:,end].*100,label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,gap_Stochastic_gen[:,end],label = "Sto-Gen", linewidth =linew)
#plot!(scenarios,gap_Stochastic_boot[:,end].*100,label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.gap .*100,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"Gap_BeforeRealization_zoomed.png")

##########################
# Display ballast water from different models
p = plot(scenarios, ones(sc)*Deterministic_Solution.ballast_weight, marker = :utriangle, xlabel = "Scenarios", 
ylabel = "Ballast weight (t)", label = "Deterministic", title = "Ballast weight for different models",
formatter=:plain)
annotate!(scenarios[1], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))
annotate!(scenarios[2], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))
annotate!(scenarios[3], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))
annotate!(scenarios[4], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))
annotate!(scenarios[5], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))

xtemp = [[],[],[],[]]
ytemp = [[],[],[],[]]
cargo_n = [[],[],[],[]]
for i in 1:sc
    if Stochastic_gen_fitted[1,i,end].status != "INFEASIBLE" # have a solution 
        push!(xtemp[1],scenarios[i])
        push!(ytemp[1],Stochastic_gen_fitted[1,i,end].ballast_weight)
        push!(cargo_n[1],Stochastic_gen_fitted[1,i,end].n_cargo_loaded)
    end
    if Stochastic_boot_fitted[1,i,end].status != "INFEASIBLE" # have a solution 
        push!(xtemp[2],scenarios[i])
        push!(ytemp[2],Stochastic_boot_fitted[1,i,end].ballast_weight)
        push!(cargo_n[2],Stochastic_boot_fitted[1,i,end].n_cargo_loaded)
    end
    #if EVP_gen_fitted[1,i,end].status != "INFEASIBLE" # have a solution 
    #    push!(xtemp[3],scenarios[i])
    #    push!(ytemp[3],EVP_gen_fitted[1,i,end].ballast_weight)
    #    push!(cargo_n[3],EVP_gen_fitted[1,i,end].n_cargo_loaded)
    #end
    #if EVP_boot_fitted[1,i,end].status != "INFEASIBLE" # have a solution 
    #    push!(xtemp[4],scenarios[i])
    #    push!(ytemp[4],EVP_boot_fitted[1,i,end].ballast_weight)
    #    push!(cargo_n[4],EVP_boot_fitted[1,i,end].n_cargo_loaded)
    #end
end
plot!(xtemp[1],ytemp[1], label = "Stochastic Gen", marker = :circle, markersize = 4)
plot!(xtemp[2],ytemp[2], label = "Stochastic Boot", marker = :circle, markersize = 4)
plot!(formatte=:plain)
#plot!(xtemp[3],ytemp[3], label = "EVP Gen", marker = :circle, markersize = 4)
#plot!(xtemp[4],ytemp[4], label = "EVP Boot", marker = :circle, markersize = 4)
for j in 1:length(xtemp)
    for i = 1:length(xtemp[j])
        annotate!(xtemp[j][i], ytemp[j][i] + 20, text(cargo_n[j][i],6))
    end
end
display(p)
savefig(plot_folder*"BallastWeight.png")

# Number of times problem was unfeasible
#EVP_gen_inf = []
#EVP_gen_fitted_inf = []
Stochastic_gen_inf = []
Stochastic_gen_fitted_inf = []
#EVP_boot_inf = []
#EVP_boot_fitted_inf = []
Stochastic_boot_inf = []
Stochastic_boot_fitted_inf = []
println("Checing problems for infeasibility")
for i in 1:sc
    for j in 1:n
        # print if not optimal
        #if EVP_gen[r,i,j].status != "OPTIMAL"
        #    push!(EVP_gen_inf,EVP_gen[r,i,j])
        #    println("##########################")
        #    println("EVP gen model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
        #    println("Status: ", EVP_gen[r,i,j].status)
        #end
        #if EVP_gen_fitted[r,i,j].status != "OPTIMAL"
        #    push!(EVP_gen_fitted_inf,EVP_gen_fitted[r,i,j])
        #    println("##########################")
        #    println("EVP gen fitted model jnfeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
        #    println("Status: ", EVP_gen_fitted[r,i,j].status)
        #end
        if Stochastic_gen[r,i,j].status != "OPTIMAL"
            push!(Stochastic_gen_inf,Stochastic_gen[r,i,j])
            println("##########################")
            println("Stochastic gen model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
            println("Status: ", Stochastic_gen[r,i,j].status)
        end
        if Stochastic_gen_fitted[r,i,j].status != "OPTIMAL"
            push!(Stochastic_gen_fitted_inf,Stochastic_gen_fitted[r,i,j])
            println("##########################")
            println("Stochastic gen fitted model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
            println("Status: ", Stochastic_gen_fitted[r,i,j].status)
        end
        #if EVP_boot[r,i,j].status != "OPTIMAL"
        #    push!(EVP_boot_inf,EVP_boot[r,i,j])
        #    println("##########################")
        #    println("EVP boot model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
        #    println("Status: ", EVP_boot[r,i,j].status)
        #end
        #if EVP_boot_fitted[r,i,j].status != "OPTIMAL"
        #    push!(EVP_boot_fitted_inf,EVP_boot_fitted[r,i,j])
        #    println("##########################")
        #    println("EVP boot fitted model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
        #    println("Status: ", EVP_boot_fitted[r,i,j].status)
        #end
        if Stochastic_boot[r,i,j].status != "OPTIMAL"
            push!(Stochastic_boot_inf,Stochastic_boot[r,i,j])
            println("##########################")
            println("Stochastic boot model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
            println("Status: ", Stochastic_boot[r,i,j].status)
        end
        if Stochastic_boot_fitted[r,i,j].status != "OPTIMAL"
            push!(Stochastic_boot_fitted_inf,Stochastic_boot_fitted[r,i,j])
            println("##########################")
            println("Stochastic boot fitted model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
            println("Status: ", Stochastic_boot_fitted[r,i,j].status)
        end
    end
end
println("##########################")
println("Number of different parameters: ", sc*n)
#println("Stochastic models + EVP models: ", length(EVP_gen)+length(EVP_boot)+length(Stochastic_gen)+length(Stochastic_boot))
println("##########################")
# EVP models - normally none of them should be infeasible
#=
for i in 1:length(EVP_gen_inf)
    println("EVP-gen Status: ", EVP_gen_inf[i].status)
    println("Index: ", i)
end
println("##########################")
for i in 1:length(EVP_boot_inf)
    println("EVP-boot Status: ", EVP_boot_inf[i].status)
    println("Index: ", i)
end
=#
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
#=
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
=#
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
println("########################")
slack_variables_boot = []
slack_variables_boot_index = []
for i in 1:sc
    foldername = "Stochastic_Bootstrap1_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        push!(slack_variables_boot_index,i)
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
println("Model indecies infeasible for Stochastic Boot: ", slack_variables_boot_index)
println("Slack variables for Stochastic Boot")
for i in 1:length(slack_variables_boot)
    for j in 1:length(slack_variables_boot[i])
        if any(slack_variables_boot[i][j] .!= 0.0)
            if j == 1
                println("In model $(slack_variables_boot_index[i]), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_boot[i][j])
            else
                println("Model $(slack_variables_boot_index[i]) and slack variable index $(j) are not zero")
                println(slack_variables_boot[i][j])
            end
        end
    end
end

# EVP Boot
#=
slack_variables_EVP_boot = []
slack_variables_EVP_boot_index = []
for i in 1:sc
    foldername = "EVP_Bootstrap1_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        push!(slack_variables_EVP_boot_index,i)
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
println("Model indecies infeasible for EVP Boot: ", slack_variables_EVP_boot_index)
println("Slack variables for EVP Boot")
for i in 1:length(slack_variables_EVP_boot)
    for j in 1:length(slack_variables_EVP_boot[i])
        if any(slack_variables_EVP_boot[i][j] .!= 0.0)
            if j == 1
                println("In model $(slack_variables_EVP_boot_index[i]), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_EVP_boot[i][j])
            else
                println("Model $(slack_variables_EVP_boot_index[i]) and slack variable index $(j) are not zero")
                println(slack_variables_EVP_boot[i][j])
            end
        end
    end
end
=#
# Stochastic Gen 
slack_variables_gen = []
slack_variables_gen_index = []
for i in 1:sc
    foldername = "Stochastic_random_sampling_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    if (temp.status == "INFEASIBLE") # solve slack model
        push!(slack_variables_gen_index,i)
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
println("Model indecies infeasible for Stochastic Gen: ", slack_variables_gen_index)
println("Slack variables for Stochastic Gen")
for i in 1:length(slack_variables_gen)
    for j in 1:length(slack_variables_gen[i])
        if any(slack_variables_gen[i][j] .!= 0.0)
            if j == 1
                println("In model $(slack_variables_gen_index[i]), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_gen[i][j])
            else
                println("Model $(slack_variables_gen_index[i]) and slack variable index $(j) are not zero")
                println(slack_variables_gen[i][j])
            end
        end
    end
end
# EVP Gen 
#=
slack_variables_EVP_gen = []
slack_variables_EVP_gen_index = []
for i in 1:sc
    foldername = "EVP_random_sampling_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        push!(slack_variables_EVP_gen_index,i)
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
            error("Model infeasible - do something to slack model")
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
println("Model indecies infeasible for EVP Gen: ", slack_variables_EVP_gen_index)
println("Slack variables for EVP Gen")
for i in 1:length(slack_variables_EVP_gen)
    for j in 1:length(slack_variables_EVP_gen[i])
        if any(abs(slack_variables_EVP_gen[i][j]) .!= 0.0)
            if j == 1
                println("In model $(slack_variables_EVP_gen_index[i]), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_EVP_gen[i][j])
            else
                println("Model $(slack_variables_EVP_gen_index[i]) and slack variable index $(j) are not zero")
                println(slack_variables_EVP_gen[i][j])
            end
        end
    end
end
=#


# Looking into limiting factor - Deck 1
println("##########################")
n_decks = length(Deterministic_problem.vessel.decks)
deck_limits = [Deterministic_problem.vessel.decks[i].weight_limit for i in 1:n_decks]
println("Deterministic Solution weight on decks:")
#secu_var = [collect(filter(x -> x.QuantileNumber == i, df_secu_quan).Variance) for i in 1:n_quantiles]
Det_sol_PlacementDecks = [filter(x -> x.deck == i, Deterministic_Solution.cargo) for i in 1:n_decks]
for i in 1:n_decks
    println("Weight Limit on deck $(i): ", deck_limits[i])
    println("Weight on deck $(i): ", sum(Det_sol_PlacementDecks[i][j].weight for j in 1:length(Det_sol_PlacementDecks[i])))
end

xtemp = [[],[],[],[]]
deck_diff = [[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]]
for i in 1:repetitions
    for j in 1:sc
        for l in 1:n
            #if EVP_gen_fitted[i,j,l].status == "TIME_LIMIT" || EVP_gen_fitted[i,j,l].status == "OPTIMAL"
            #    EVP_gen_fitted_PlacementDecks = [filter(x -> x.deck == k, EVP_gen_fitted[i,j,l].cargo) for k in 1:n_decks]
            #    println("#################")
            #    push!(xtemp[1],scenarios[j])
            #    for k in 1:n_decks
            #        push!(deck_diff[k][1], deck_limits[k]-sum(EVP_gen_fitted_PlacementDecks[k][m].weight for m in 1:length(EVP_gen_fitted_PlacementDecks[k])))
            #        println("EVP gen fitted model $(i), scenarios $(scenarios[j]), unknown cargo $(n_unknown[l]) weight-difference on deck $(k): ", deck_limits[k]-sum(EVP_gen_fitted_PlacementDecks[k][m].weight for m in 1:length(EVP_gen_fitted_PlacementDecks[k])))
            #    end
            #else # not feasible, probably due to deck 1
            #    push!(xtemp[1],scenarios[j])
            #    EVP_gen_PlacementDecks = [filter(x -> x.deck == k,EVP_gen[i,j,l].cargo) for k in 1:n_decks]
                # calculate how much over limit - same as slack variables?
            #    for k in 1:n_decks
            #        actual_weight = [filter(x -> x.id == EVP_gen_PlacementDecks[k][i].id, Deterministic_problem.cargo.items)[1].weight for i in 1:length(EVP_gen_PlacementDecks[k])]
            #        push!(deck_diff[k][1], deck_limits[k]-sum(actual_weight))
            #    end
            #end
            #=
            if EVP_boot_fitted[i,j,l].status == "TIME_LIMIT" || EVP_boot_fitted[i,j,l].status == "OPTIMAL"
                EVP_boot_fitted_PlacementDecks = [filter(x -> x.deck == k, EVP_boot_fitted[i,j,l].cargo) for k in 1:n_decks]
                println("#################")
                push!(xtemp[2],scenarios[j])
                for k in 1:n_decks
                    push!(deck_diff[k][2], deck_limits[k]-sum(EVP_boot_fitted_PlacementDecks[k][m].weight for m in 1:length(EVP_boot_fitted_PlacementDecks[k])))
                    println("EVP boot fitted model $(i), scenarios $(scenarios[j]), unknown cargo $(n_unknown[l]) weight on deck $(k): ", deck_limits[k]-sum(EVP_boot_fitted_PlacementDecks[k][m].weight for m in 1:length(EVP_boot_fitted_PlacementDecks[k])))
                end
            else # not feasible, probably due to deck 1
                push!(xtemp[2],scenarios[j])
                EVP_boot_PlacementDecks = [filter(x -> x.deck == k,EVP_boot[i,j,l].cargo) for k in 1:n_decks]
                # calculate how much over limit - same as slack variables?
                for k in 1:n_decks
                    actual_weight = [filter(x -> x.id == EVP_boot_PlacementDecks[k][i].id, Deterministic_problem.cargo.items)[1].weight for i in 1:length(EVP_boot_PlacementDecks[k])]
                    push!(deck_diff[k][2], deck_limits[k]-sum(actual_weight))
                end
            end
            =#
            if Stochastic_gen_fitted[i,j,l].status == "TIME_LIMIT" || Stochastic_gen_fitted[i,j,l].status == "OPTIMAL"
                Stochastic_gen_fitted_PlacementDecks = [filter(x -> x.deck == k, Stochastic_gen_fitted[i,j,l].cargo) for k in 1:n_decks]
                println("#################")
                push!(xtemp[3],scenarios[j])
                for k in 1:n_decks
                    push!(deck_diff[k][3], deck_limits[k]-sum(Stochastic_gen_fitted_PlacementDecks[k][m].weight for m in 1:length(Stochastic_gen_fitted_PlacementDecks[k])))
                    println("Stochastic gen fitted model $(i), scenarios $(scenarios[j]), unknown cargo $(n_unknown[l]) weight on deck $(k): ", deck_limits[k]-sum(Stochastic_gen_fitted_PlacementDecks[k][m].weight for m in 1:length(Stochastic_gen_fitted_PlacementDecks[k])))
                end
            else # not feasible, probably due to deck 1
                push!(xtemp[3],scenarios[j])
                Stochastic_gen_PlacementDecks = [filter(x -> x.deck == k,Stochastic_gen[i,j,l].cargo) for k in 1:n_decks]
                # calculate how much over limit - same as slack variables?
                for k in 1:n_decks
                    actual_weight = [filter(x -> x.id == Stochastic_gen_PlacementDecks[k][i].id, Deterministic_problem.cargo.items)[1].weight for i in 1:length(Stochastic_gen_PlacementDecks[k])]
                    push!(deck_diff[k][3], deck_limits[k]-sum(actual_weight))
                end
            end
            if Stochastic_boot_fitted[i,j,l].status == "TIME_LIMIT" || Stochastic_boot_fitted[i,j,l].status == "OPTIMAL"
                Stochastic_boot_fitted_PlacementDecks = [filter(x -> x.deck == k, Stochastic_boot_fitted[i,j,l].cargo) for k in 1:n_decks]
                println("#################")
                push!(xtemp[4],scenarios[j])
                for k in 1:n_decks
                    push!(deck_diff[k][4], deck_limits[k]-sum(Stochastic_boot_fitted_PlacementDecks[k][m].weight for m in 1:length(Stochastic_boot_fitted_PlacementDecks[k])))
                    println("Stochastic boot fitted model $(i), scenarios $(scenarios[j]), unknown cargo $(n_unknown[l]) weight on deck $(k): ", deck_limits[k]-sum(Stochastic_boot_fitted_PlacementDecks[k][m].weight for m in 1:length(Stochastic_boot_fitted_PlacementDecks[k])))
                end
            else # not feasible, probably due to deck 1
                push!(xtemp[4],scenarios[j])
                Stochastic_boot_PlacementDecks = [filter(x -> x.deck == k,Stochastic_boot[i,j,l].cargo) for k in 1:n_decks]
                # calculate how much over limit - same as slack variables?
                for k in 1:n_decks
                    actual_weight = [filter(x -> x.id == Stochastic_boot_PlacementDecks[k][i].id, Deterministic_problem.cargo.items)[1].weight for i in 1:length(Stochastic_boot_PlacementDecks[k])]
                    push!(deck_diff[k][4], deck_limits[k]-sum(actual_weight))
                end
            end
        end
    end
end
# plot for decks
p = []
for i in 1:n_decks
    pl = plot(xtemp[3],deck_diff[i][3], label = "Sto-gen", marker = :circle, markersize = 4,
    xlabel="Scenarios", ylabel="Weight diff. (t)",
    title = "Unused Weight Capacity on Deck $(i)\n for different models")
    #plot!(xtemp[2],deck_diff[i][2], label = "EVP-boot", marker = :circle, markersize = 4)
    #plot!(xtemp[3],deck_diff[i][3], label = "Sto-gen", marker = :circle, markersize = 4)
    plot!(xtemp[4],deck_diff[i][4], label = "Sto-boot", marker = :circle, markersize = 4)
    savefig(plot_folder*"Deck$(i)_WeightDiff.png")
    push!(p,pl)
end
for i in p
    display(i)
end
