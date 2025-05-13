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
using StatsPlots

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
include("src/representation/VarianceOfWeight.jl")


# load data from problems
# Change to get different problems
#HPC_folder = "Finlandia_01_04_09_41_38"
#HPC_folder = "Finlandia_07_04_15"
#plot_folder = "Plots/Data/07_04_15/"
HPC_folder = "Finlandia_09_04_09"
plot_folder = "Plots/Data/Finlandia_09_04_09/"

# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end

repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
println("Extra info about test: ", note)
sc = length(scenarios)
n = length(n_unknown)
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous"
Deterministic_problem = load_data(problemname1,problemname2,problemname3)
EVP_problem_gen = Array{Any}(undef, repetitions, sc,n)
Stochastic_problem_gen = Array{Any}(undef, repetitions, sc,n)
EVP_problem_boot = Array{Any}(undef, repetitions, sc,n)
Stochastic_problem_boot = Array{Any}(undef, repetitions, sc,n)

for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # Load EVP
            foldername = "EVP_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Problem"
            EVP_problem_gen[i,j,k] = get_deterministic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
            foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Problem"
            EVP_problem_boot[i,j,k] = get_deterministic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
            # Load Stochastic
            foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem_gen[i,j,k] = get_stochastic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
            foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem_boot[i,j,k] = get_stochastic_problem(foldername,
            filename,HPC_folder,problemname1,problemname2,problemname3)
        end
    end
end
#################################
# Plot Cargo weights for some scenarios for a problem 
plots_boot = plot_cargo_weights(Stochastic_problem_boot[1,1,1],[i for i =1:scenarios[1]]) # plot for multiple scenarios
plots_gen = plot_cargo_weights(Stochastic_problem_gen[1,1,1],[i for i =1:scenarios[1]])
xplots = 2
yplots = 2
plot(plots_boot[1],plots_boot[2],plots_boot[3],plots_boot[4], layout=(xplots,yplots)) 
savefig(plot_folder*"Cargo_distribution_boot_scenario_1_4.png")
plot(plots_boot[5],plots_boot[6],plots_boot[7],plots_boot[8], layout=(xplots,yplots)) 
savefig(plot_folder*"Cargo_distribution_boot_scenario_5_8.png")
plot(plots_boot[9],plots_boot[10], layout=(1,2)) 
savefig(plot_folder*"Cargo_distribution_boot_scenario_9_10.png")
plot(plots_gen[1],plots_gen[2],plots_gen[3],plots_gen[4], layout=(xplots,yplots))
savefig(plot_folder*"Cargo_distribution_gen_scenario_1_4.png")
plot(plots_gen[5],plots_gen[6],plots_gen[7],plots_gen[8], layout=(xplots,yplots)) 
savefig(plot_folder*"Cargo_distribution_gen_scenario_5_8.png")
plot(plots_gen[9],plots_gen[10], layout=(1,2))
savefig(plot_folder*"Cargo_distribution_gen_scenario_9_10.png")
 # plot for original problem
OG_plot = plot_cargo_OG(Deterministic_problem)
EVP_plot_boot = plot_cargo_OG(EVP_problem_boot[1,1,end],false) # plot for EVP
#plot!(EVP_plot_boot[1], title = "Cargo weight distribution for EVP - Boot")
#plot!(EVP_plot_boot[2], title = "Cargo weight distribution for EVP - Boot")
EVP_plot_gen = plot_cargo_OG(EVP_problem_gen[1,1,end],false)
#plot!(EVP_plot_gen[1], title = "Cargo weight distribution for EVP - Gen")
#plot!(EVP_plot_gen[2], title = "Cargo weight distribution for EVP - Gen")

display(OG_plot[1])
display(OG_plot[2])
display(EVP_plot_boot[1])
display(EVP_plot_boot[2])
display(EVP_plot_gen[1])
display(EVP_plot_gen[2])
# Save plots
savefig(OG_plot[1],plot_folder*"Cargo_distribution_Determinstic_1.png")
savefig(OG_plot[2],plot_folder*"Cargo_distribution_Determinstic_2.png")
savefig(EVP_plot_boot[1],plot_folder*"Cargo_distribution_EVP_boot_1.png")
savefig(EVP_plot_boot[2],plot_folder*"Cargo_distribution_EVP_boot_2.png")
savefig(EVP_plot_gen[1],plot_folder*"Cargo_distribution_EVP_gen_1.png")
savefig(EVP_plot_gen[2],plot_folder*"Cargo_distribution_EVP_gen_2.png")
# total weight of cargo
println("Total weight of cargo: ", sum([cargo.weight for cargo in Deterministic_problem.cargo]))
println("Average total weight of cargo, stochastic problem - Gen: ", sum(sum([cargo.weight for cargo in Stochastic_problem_boot[1,1,end].cargo.items[i]]) for i in 1:scenarios[1])/scenarios[1])
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            println("####################")
            for l in 1:scenarios[j]
                println("Total weight of cargo - boot, repetition: ", i, ", number of scenarios: ", scenarios[j],
                ", number of unknown weights: ", n_unknown[k], ": ",
                sum([cargo.weight for cargo in Stochastic_problem_boot[i,j,k].cargo.items[l]]))
            end
        end
    end
end

################################
# Plot for total weight in each scenario - Gen 
# Change this if want problems with a different number of scenarios
for l in 1:sc
    n_sc = l
    total_weights_sto = Array{Any}(undef, scenarios[n_sc],n)
    for i in 1:n
        sto_pro = Stochastic_problem_gen[1,n_sc,i]
        for j in 1:scenarios[n_sc]
            total_weights_sto[j,i] = sum([cargo.weight for cargo in sto_pro.cargo.items[j]])
        end
    end
    local total_weight_det = ones(Stochastic_problem_gen[1,n_sc,end].scenarios)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
    local p = plot(total_weights_sto[:,1],xlabel = "Scenarios", ylabel = "Total cargo weight (t)", 
    label = "Sto_pro, n: $(n_unknown[1])",
    title = "Total weight for the stochastic problem - Gen,\nscenarios: $(scenarios[n_sc])")
    for i in 2:n
        plot!(total_weights_sto[:,i], label = "Sto_pro, n: $(n_unknown[i])")
    end
    plot!(total_weight_det, label = "Det_pro")
    display(p)
    #savefig(p, plot_folder*"Total_weight_stochastic_problem_gen_scenario_$(scenarios[n_sc]).png")
end
# NB! 
# Clearly shows that our generation of weight is generally 
# higher than the deterministic problem

# Display EVP data
    #n_sc = l
    total_weights_EVP = []
    for i in 1:sc
        push!(total_weights_EVP, sum([cargo.weight for cargo in EVP_problem_gen[1,i,end].cargo]))
    end
    total_weight_det = ones(sc)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
    p = plot(scenarios,total_weights_EVP,xlabel = "Scenarios", ylabel = "Total cargo weight (t)", label = "EVP", 
    title = "Total weight for the EVP - Gen.")
    plot!(scenarios,total_weight_det, label = "Det_pro")
    display(p)
    savefig(p, plot_folder*"Total_weight_EVP_gen.png")
# Clearly shows that our generation of weight is generally 
# higher than the deterministic problem

# Plot for total weight in each scenario - Boot
for l in 1:sc
    n_sc = l
    total_weights_sto = Array{Any}(undef, scenarios[n_sc],n)
    for i in 1:n
        sto_pro = Stochastic_problem_boot[1,n_sc,i]
        for j in 1:scenarios[n_sc]
            total_weights_sto[j,i] = sum([cargo.weight for cargo in sto_pro.cargo.items[j]])
        end
    end
    local total_weight_det = ones(Stochastic_problem_boot[1,n_sc,end].scenarios)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
    local p = plot(total_weights_sto[:,1],xlabel = "Scenarios", ylabel = "Total cargo weight (t)", 
    label = "Sto_pro, n: $(n_unknown[1])",
    title = "Total weight for the stochastic problem - Boot.")
    for i in 2:n
        plot!(total_weights_sto[:,i], label = "Sto_pro, n: $(n_unknown[i])")
    end
    plot!(total_weight_det, label = "Det_pro")
    display(p)
    #savefig(p, plot_folder*"Total_weight_stochastic_problem_boot_scenario_$(scenarios[n_sc]).png")
end
# NB! 
# Shows Boot-scenario generation generally is lower than the actual weight

# Display EVP data
total_weights_EVP = []
for i in 1:sc
    push!(total_weights_EVP, sum([cargo.weight for cargo in EVP_problem_boot[1,i,end].cargo]))
end
total_weight_det = ones(sc)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
p = plot(scenarios,total_weights_EVP,xlabel = "scenarios", ylabel = "Total cargo weight (t)", label = "EVP", 
title = "Total weight for the EVP - Boot.")
plot!(scenarios,total_weight_det, label = "Det_pro")
savefig(p, plot_folder*"Total_weight_EVP_boot.png")
display(p)
# Clearly shows that our generation of weight is generally 
# higher than the deterministic problem


# Other methods for generating scenarios...


########################################
# Plot Historic data
plot_folder_historic = "Plots/Data/Historic/"
# Create folder for plots
if !isdir(plot_folder_historic)
    mkpath(plot_folder_historic)
end
# load historic data 
n_quantiles = 4
cur_path = @__DIR__
file_path = joinpath(cur_path,"data","CargoWeights.csv")
df = load_Weight_Variance_data(file_path)
df_secu, df_trailer = weight_difference_info(df,false)
remove_outliers = false
df_secu_quan, = seperate_data_into_quantiles(df_secu,n_quantiles,remove_outliers)
df_trailer_quan, = seperate_data_into_quantiles(df_trailer,n_quantiles,remove_outliers)
normal_var_secu = (df_secu_quan.Variance.-mean(df_secu_quan.Variance))./std(df_secu_quan.Variance)
qqplot(normal_var_secu, Normal(0, 1), title="Trucks Q-Q Plot", xlabel="Theoretical Quantiles", ylabel="Truck Quantiles")
savefig(plot_folder_historic*"QQPlot_secu.png")
normal_var_truck = (df_trailer_quan.Variance.-mean(df_trailer_quan.Variance))./std(df_trailer_quan.Variance)
qqplot(normal_var_truck, Normal(0, 1), title="Trucks Q-Q Plot", xlabel="Theoretical Quantiles", ylabel="Truck Quantiles")
savefig(plot_folder_historic*"QQPlot_truck.png")
# Removing outliers
remove_outliers = true
df_secu_quan, = seperate_data_into_quantiles(df_secu,n_quantiles,remove_outliers)
df_trailer_quan, = seperate_data_into_quantiles(df_trailer,n_quantiles,remove_outliers)

# Booked weight
nbins = 20
test = fit(Histogram, df_secu_quan.CountBookedWeight, nbins=nbins)
binedges = test.edges[1]
histogram(df_secu_quan.CountBookedWeight,bins = nbins, xlabel="Weight", 
ylabel="Frequency", title="Secu-boxes booked weight distribution", legend=false)
xticks!(binedges)
plot!(xtickfontsize=6,formatter=:plain)
savefig(plot_folder_historic*"Weight_booked_Secu.png")
# Er den her bedre?
histogram(df_secu_quan.CountBookedWeight, xlabel="Weight", 
ylabel="Frequency", title="SECU booked weight distribution", legend=false)
plot!(formatter=:plain)

nbins = 30
test = fit(Histogram, df_trailer_quan.CountBookedWeight, nbins=nbins)
binedges = test.edges[1]
histogram(df_trailer_quan.CountBookedWeight, bins=nbins, xlabel="Weight",
ylabel="Frequency", title="Trailer booked weight distribution", legend=false)
xticks!(binedges)
plot!(xtickfontsize=4,formatter=:plain)
savefig(plot_folder_historic*"Weight_booked_trailer.png")
# Er den her bedre?
histogram(df_trailer_quan.CountBookedWeight, xlabel="Weight",
ylabel="Frequency", title="Trailer booked weight distribution", legend=false)
plot!(formatter=:plain)

# Weight variance
nbins = 20
test = fit(Histogram, df_secu_quan.Variance, nbins=nbins)
binedges = test.edges[1]
histogram(df_secu_quan.Variance,bins = nbins, xlabel="Weight (kg)", 
ylabel="Frequency", title="Secu-boxes weight difference distribution", legend=false)
xticks!(binedges)
plot!(xtickfontsize=5,formatter=:plain)
#savefig(plot_folder_historic*"Weight_variance_Secu.png")
# Er det bedre ikke at bestemme bin-width?
histogram(df_secu_quan.Variance./1000, xlabel="Weight (t.)", 
ylabel="Frequency", title="Secu-boxes weight difference distribution", legend=false)
plot!(xtickfontsize=6,formatter=:plain)
savefig(plot_folder_historic*"Weight_variance_Secu.png")

nbins = 20
test = fit(Histogram, df_trailer_quan.Variance, nbins=nbins)
binedges = test.edges[1]
histogram(df_trailer_quan.Variance, bins=nbins, xlabel="Weight",
ylabel="Frequency", title="Trucks weight difference distribution", legend=false)
xticks!(binedges)
plot!(xtickfontsize=4,formatter=:plain)
#savefig(plot_folder_historic*"Weight_variance_trailer.png")
# Bedre?
histogram(df_trailer_quan.Variance./1000, bins=nbins, xlabel="Weight (t.)",
ylabel="Frequency", title="Trucks weight difference distribution", legend=false)
plot!(xtickfontsize=6,formatter=:plain)
savefig(plot_folder_historic*"Weight_variance_trailer.png")

secu_var = [collect(filter(x -> x.QuantileNumber == i, df_secu_quan).Variance) for i in 1:n_quantiles]
trailer_var = [collect(filter(x -> x.QuantileNumber == i, df_trailer_quan).Variance) for i in 1:n_quantiles]
p = []
for i in 1:n_quantiles
    push!(p, histogram(secu_var[i], bins=30, xlabel="Weight variance", 
    ylabel="Frequency", title="SECU weight variance, quantile: $i",titlefont=font(10),tickfontsize=5, legend=false, formatter=:plain))
end
xplots = Int(ceil(n_quantiles/2))
yplots = Int(ceil(n_quantiles/xplots))
plot(p..., layout=(xplots,yplots)) # Good plot I think
savefig(plot_folder_historic*"Weight_variance_secu_quantiles.png")
p = []
for i in 1:n_quantiles
    push!(p, histogram(trailer_var[i], bins=30, xlabel="Weight variance", 
    ylabel="Frequency", title="Trailer weight variance, quantile: $i",titlefont=font(10),tickfontsize=5, legend=false, formatter=:plain))
end
plot(p..., layout=(xplots,yplots)) # Good plot I think
savefig(plot_folder_historic*"Weight_variance_trailer_quantiles.png")



q = [1,2,3,4]#,5,6,7,8]
p = []
cargoC_boot = Stochastic_problem_boot[1,1,1].cargo.items[1]
cargoC_gen = Stochastic_problem_gen[1,1,1].cargo.items[1]
for i in 1:length(q)
    df_secu_quan_temp, = seperate_data_into_quantiles(df_secu,q[i],true)
    df_trailer_quan_temp, = seperate_data_into_quantiles(df_trailer,q[i],true)
    test_plot = scenario_distribution_with_normal(cargoC_boot,df_secu_quan_temp,df_trailer_quan_temp,q[i])
    push!(p,test_plot)
    display(test_plot[1])
    display(test_plot[2])
    #println("i: $(i), secu: ", secu_weight)
    #println("i: $(i), trailer: ", trailer_weight)
end

