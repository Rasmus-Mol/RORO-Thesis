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
include("src/representation/VarianceOfWeight.jl")


# load data from problems
# Chance to get different problems
HPC_folder = "Finlandia_01_04_09_41_38"

repetitions, scenarios, n_unknown, time_limit = get_HPC_data(HPC_folder)
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
plots_gen = plot_cargo_weights(Stochastic_problem_boot[1,1,end],[1,2,3]) # plot for multiple scenarios
for i in 1:length(plots_gen)
    display(plots_gen[i])
end
OG_plot = plot_cargo_OG(Deterministic_problem) # plot for original problem
EVP_plot = plot_cargo_OG(EVP_problem_boot[1,1,end],false) # plot for EVP
display(OG_plot[1])
display(OG_plot[2])
display(EVP_plot[1])
display(EVP_plot[2])
# total weight of cargo
println("Total weight of cargo: ", sum([cargo.weight for cargo in Deterministic_problem.cargo]))
println("Average total weight of cargo, stochastic problem - Gen: ", sum(sum([cargo.weight for cargo in Stochastic_problem_boot[1,1,end].cargo.items[i]]) for i in 1:scenarios[1])/scenarios[1])
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            println("####################")
            for l in 1:scenarios[i]
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
n_sc = 1
total_weights_sto = Array{Any}(undef, scenarios[n_sc],n)
for i in 1:n
    sto_pro = Stochastic_problem_gen[1,n_sc,i]
    for j in 1:scenarios[n_sc]
        total_weights_sto[j,i] = sum([cargo.weight for cargo in sto_pro.cargo.items[j]])
    end
end
total_weight_det = ones(Stochastic_problem_gen[1,n_sc,end].scenarios)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
p = plot(total_weights_sto[:,1],xlabel = "Scenario", ylabel = "Total cargo weight (t)", 
label = "Sto_pro, n: $(n_unknown[1])",
title = "Total weight for the stochastic problem - Gen,\nscenarios: $(scenarios[n_sc])")
for i in 2:n
    plot!(total_weights_sto[:,i], label = "Sto_pro, n: $(n_unknown[i])")
end
plot!(total_weight_det, label = "Det_pro")
# NB! 
# Clearly shows that our generation of weight is generally 
# higher than the deterministic problem
display(p)

# Display EVP data
n_sc = 1
total_weights_EVP = []
for i in 1:n
    push!(total_weights_EVP, sum([cargo.weight for cargo in EVP_problem_gen[1,n_sc,i].cargo]))
end
total_weight_det = ones(n)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
p = plot(n_unknown,total_weights_EVP,xlabel = "n_unknown", ylabel = "Total cargo weight (t)", label = "EVP", 
title = "Total weight for the EVP - Boot,\n scenarios: $(scenarios[n_sc]), ")
plot!(n_unknown,total_weight_det, label = "Det_pro")
# Clearly shows that our generation of weight is generally 
# higher than the deterministic problem
display(p)

# Plot for total weight in each scenario - Boot
# Change this if want problems with a different number of scenarios
n_sc = 1
total_weights_sto = Array{Any}(undef, scenarios[n_sc],n)
for i in 1:n
    sto_pro = Stochastic_problem_boot[1,n_sc,i]
    for j in 1:scenarios[n_sc]
        total_weights_sto[j,i] = sum([cargo.weight for cargo in sto_pro.cargo.items[j]])
    end
end
total_weight_det = ones(Stochastic_problem_boot[1,n_sc,end].scenarios)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
p = plot(total_weights_sto[:,1],xlabel = "Scenario", ylabel = "Total cargo weight (t)", 
label = "Sto_pro, n: $(n_unknown[1])",
title = "Total weight for the stochastic problem - Boot,\nscenarios: $(scenarios[n_sc])")
for i in 2:n
    plot!(total_weights_sto[:,i], label = "Sto_pro, n: $(n_unknown[i])")
end
plot!(total_weight_det, label = "Det_pro")
# NB! 
# Shows Boot-scenario generation generally is lower than the actual weight
display(p)

# Display EVP data
n_sc = 1
total_weights_EVP = []
for i in 1:n
    push!(total_weights_EVP, sum([cargo.weight for cargo in EVP_problem_boot[1,n_sc,i].cargo]))
end
total_weight_det = ones(n)*sum([cargo.weight for cargo in Deterministic_problem.cargo])
p = plot(n_unknown,total_weights_EVP,xlabel = "n_unknown", ylabel = "Total cargo weight (t)", label = "EVP", 
title = "Total weight for the EVP - Boot,\n scenarios: $(scenarios[n_sc]), ")
plot!(n_unknown,total_weight_det, label = "Det_pro")
# Clearly shows that our generation of weight is generally 
# higher than the deterministic problem
display(p)

# Other methods for generating scenaros...


########################################
# Plot Historic data
# load historic data 
n_quantiles = 4
cur_path = @__DIR__
file_path = joinpath(cur_path,"data","CargoWeights.csv")
df = load_Weight_Variance_data(file_path)
df_secu, df_trailer = weight_difference_info(df,false)
df_secu_quan, = seperate_data_into_quantiles(df_secu,n_quantiles,false)
df_trailer_quan, = seperate_data_into_quantiles(df_trailer,n_quantiles,false)

# Booked weight
nbins = 20
test = fit(Histogram, df_secu_quan.CountBookedWeight, nbins=nbins)
binedges = test.edges[1]
histogram(df_secu_quan.CountBookedWeight,bins = nbins, xlabel="Weight", 
ylabel="Frequency", title="SECU booked weight distribution", legend=false)
xticks!(binedges)
plot!(xtickfontsize=6,formatter=:plain)
savefig("Plots/Data/Weight_booked_Secu.png")
nbins = 30
test = fit(Histogram, df_trailer_quan.CountBookedWeight, nbins=nbins)
binedges = test.edges[1]
histogram(df_trailer_quan.CountBookedWeight, bins=nbins, xlabel="Weight",
ylabel="Frequency", title="Trailer booked weight distribution", legend=false)
xticks!(binedges)
plot!(xtickfontsize=4,formatter=:plain)
savefig("Plots/Data/Weight_booked_trailer.png")
# Weight variance
nbins = 20
test = fit(Histogram, df_secu_quan.CountBookedWeight, nbins=nbins)
binedges = test.edges[1]
histogram(df_secu_quan.Variance,bins = nbins, xlabel="Weight", 
ylabel="Frequency", title="SECU variance distribution", legend=false)
xticks!(binedges)
plot!(xtickfontsize=6,formatter=:plain)
savefig("Plots/Data/Weight_variance_Secu.png")
nbins = 30
test = fit(Histogram, df_trailer_quan.CountBookedWeight, nbins=nbins)
binedges = test.edges[1]
histogram(df_trailer_quan.Variance, bins=nbins, xlabel="Weight",
ylabel="Frequency", title="Trailer variance distribution", legend=false)
xticks!(binedges)
plot!(xtickfontsize=4,formatter=:plain)
savefig("Plots/Data/Weight_variance_trailer.png")

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
savefig("Plots/Data/Weight_variance_secu_quantiles.png")
p = []
for i in 1:n_quantiles
    push!(p, histogram(trailer_var[i], bins=30, xlabel="Weight variance", 
    ylabel="Frequency", title="Trailer weight variance, quantile: $i",titlefont=font(10),tickfontsize=5, legend=false, formatter=:plain))
end
plot(p..., layout=(xplots,yplots)) # Good plot I think
savefig("Plots/Data/Weight_variance_trailer_quantiles.png")
