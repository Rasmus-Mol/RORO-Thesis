include("packages_and_files.jl")

# New analysis
plot_folder = "Plots/Data/Scenarios_new/"
# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end
HPC_folder = "Finlandia_no_cars_heavy_100_Scenario_test_20_05_18"
repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
test_problem = Finlandia_test[8]
sc = length(scenarios)
n = length(n_unknown)
scenarios_test = [100, 200, 400] # Hardcoded in HPC-scripts
no_test = length(scenarios_test)
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", test_problem, "hazardous"
Deterministic_problem = load_data(problemname1, problemname2, problemname3)
Stochastic_problem_gen = Array{Any}(undef, repetitions, sc, n)
Stochastic_problem_boot = Array{Any}(undef, repetitions, sc, n)
Stochastic_problem_gen_Naive = Array{Any}(undef, repetitions, sc, n, no_test)
Stochastic_problem_boot_Naive = Array{Any}(undef, repetitions, sc, n, no_test)
Stochastic_problem_gen_Heuristic = Array{Any}(undef, repetitions, sc, n, no_test)
Stochastic_problem_boot_Heuristic = Array{Any}(undef, repetitions, sc, n, no_test)

Stochastic_solution_gen = Array{Any}(undef, repetitions, sc, n)
Stochastic_solution_boot = Array{Any}(undef, repetitions, sc, n)
Stochastic_solution_gen_Naive = Array{Any}(undef, repetitions, sc, n, no_test)
Stochastic_solution_boot_Naive = Array{Any}(undef, repetitions, sc, n, no_test)
Stochastic_solution_gen_Heuristic = Array{Any}(undef, repetitions, sc, n, no_test)
Stochastic_solution_boot_Heuristic = Array{Any}(undef, repetitions, sc, n, no_test)

for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # Load Stochastic - no reduction
            foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem_gen[i, j, k] = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, problemname2, problemname3)
            filename = "Stochastic_Solution"
            Stochastic_solution_gen[i, j, k] = get_solution_stochastic(foldername,
                filename, HPC_folder)
            foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem_boot[i, j, k] = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, problemname2, problemname3)
            filename = "Stochastic_Solution"
            Stochastic_solution_boot[i, j, k] = get_solution_stochastic(foldername,
                filename, HPC_folder)
            for l in 1:no_test
                # Load Stochastic - Naive
                foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)_Naive_sc$(scenarios_test[l])"
                filename = "Stochastic_Problem"
                Stochastic_problem_gen_Naive[i, j, k, l] = get_stochastic_problem(foldername,
                    filename, HPC_folder, problemname1, problemname2, problemname3)
                filename = "Stochastic_Solution"
                Stochastic_solution_gen_Naive[i, j, k, l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
                filename = "Stochastic_Problem"
                foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)_Naive_sc$(scenarios_test[l])"
                Stochastic_problem_boot_Naive[i, j, k, l] = get_stochastic_problem(foldername,
                    filename, HPC_folder, problemname1, problemname2, problemname3)
                filename = "Stochastic_Solution"
                Stochastic_solution_boot_Naive[i, j, k, l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)

                # Load Stochastic - Heuristic
                foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)_Heuristic_sc$(scenarios_test[l])"
                filename = "Stochastic_Problem"
                Stochastic_problem_gen_Heuristic[i, j, k, l] = get_stochastic_problem(foldername,
                    filename, HPC_folder, problemname1, problemname2, problemname3)
                filename = "Stochastic_Solution"
                Stochastic_solution_gen_Heuristic[i, j, k, l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
                foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)_Heuristic_sc$(scenarios_test[l])"
                filename = "Stochastic_Problem"
                Stochastic_problem_boot_Heuristic[i, j, k, l] = get_stochastic_problem(foldername,
                    filename, HPC_folder, problemname1, problemname2, problemname3)
                filename = "Stochastic_Solution"
                Stochastic_solution_boot_Heuristic[i, j, k, l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
            end
        end
    end
end
# Results
# Cargo loaded
cargo_loaded_gen = [sol.n_cargo_loaded for sol in Stochastic_solution_gen]
cargo_loaded_boot = [sol.n_cargo_loaded for sol in Stochastic_solution_boot]
cargo_loaded_gen_Naive = [sol.n_cargo_loaded for sol in Stochastic_solution_gen_Naive]
cargo_loaded_boot_Naive = [sol.n_cargo_loaded for sol in Stochastic_solution_boot_Naive]
cargo_loaded_gen_Heuristic = [sol.n_cargo_loaded for sol in Stochastic_solution_gen_Heuristic]
cargo_loaded_boot_Heuristic = [sol.n_cargo_loaded for sol in Stochastic_solution_boot_Heuristic]
# Gen
p = plot(scenarios,cargo_loaded_gen[1,:,1], label="No reduction", xlabel="Scenarios", ylabel="Cargo loaded", title="Cargo loaded for different scenarios,\n Uniform random sampling.")
#plot!(scenarios,cargo_loaded_gen[2,:,1], label="Gen rep 2")
for i in 1:no_test
    plot!(scenarios, cargo_loaded_gen_Naive[1,:,1,i], label="M 1,\n $(scenarios_test[i]) samples")
    plot!(scenarios, cargo_loaded_gen_Heuristic[1,:,1,i], label="M 2,\n $(scenarios_test[i]) samples")
end
display(p)
p = plot(scenarios,cargo_loaded_gen[2,:,1], label="No reduction",legend=:bottomleft, xlabel="Scenarios", ylabel="Cargo loaded", title="Cargo loaded for different scenarios,\n Uniform random sampling.")
#plot!(scenarios,cargo_loaded_gen[2,:,1], label="Gen rep 2")
for i in 1:no_test
    plot!(scenarios, cargo_loaded_gen_Naive[2,:,1,i], label="M 1,\n $(scenarios_test[i]) samples")
    plot!(scenarios, cargo_loaded_gen_Heuristic[2,:,1,i], label="M 2,\n $(scenarios_test[i]) samples")
end
display(p)
# Boot
p = plot(scenarios,cargo_loaded_boot[1,:,1], label="No reduction", xlabel="Scenarios", ylabel="Cargo loaded", title="Cargo loaded for different scenarios,\n sampling method 2.")
#plot!(scenarios,cargo_loaded_gen[2,:,1], label="Gen rep 2")
for i in 1:no_test
    plot!(scenarios, cargo_loaded_boot_Naive[1,:,1,i], label="M 1,\n $(scenarios_test[i]) samples")
    plot!(scenarios, cargo_loaded_boot_Heuristic[1,:,1,i], label="M 2,\n $(scenarios_test[i]) samples")
end
display(p)
p = plot(scenarios,cargo_loaded_boot[2,:,1], label="No reduction", xlabel="Scenarios", ylabel="Cargo loaded", title="Cargo loaded for different scenarios,\n sampling method 2.")
#plot!(scenarios,cargo_loaded_gen[2,:,1], label="Gen rep 2")
for i in 1:no_test
    plot!(scenarios, cargo_loaded_boot_Naive[2,:,1,i], label="M 1,\n $(scenarios_test[i]) samples")
    plot!(scenarios, cargo_loaded_boot_Heuristic[2,:,1,i], label="M 2,\n $(scenarios_test[i]) samples")
end
display(p)

# Gap
gap_gen = [sol.gap for sol in Stochastic_solution_gen]
gap_boot = [sol.gap for sol in Stochastic_solution_boot]
gap_gen_Naive = [sol.gap for sol in Stochastic_solution_gen_Naive]
gap_boot_Naive = [sol.gap for sol in Stochastic_solution_boot_Naive]
gap_gen_Heuristic = [sol.gap for sol in Stochastic_solution_gen_Heuristic]
gap_boot_Heuristic = [sol.gap for sol in Stochastic_solution_boot_Heuristic]
# Gen
println("#########################")
pretty_table(gap_gen[:,:,1],title = "Uniform random sampling - no reduction", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_gen_Naive[:,:,1,1],title = "Uniform random sampling - method 1, $(scenarios_test[1])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_gen_Naive[:,:,1,2],title = "Uniform random sampling - method 1, $(scenarios_test[2])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_gen_Naive[:,:,1,3],title = "Uniform random sampling - method 1, $(scenarios_test[3])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_gen_Heuristic[:,:,1,1],title = "Uniform random sampling - method 2, $(scenarios_test[1])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_gen_Heuristic[:,:,1,2],title = "Uniform random sampling - method 2, $(scenarios_test[2])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_gen_Heuristic[:,:,1,3],title = "Uniform random sampling - method 2, $(scenarios_test[3])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
# Boot
println("#########################")
pretty_table(gap_boot[:,:,1],title = "Sampling 2 - no reduction", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_boot_Naive[:,:,1,1],title = "Sampling 2 - method 1, $(scenarios_test[1])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_boot_Naive[:,:,1,2],title = "Sampling 2 - method 1, $(scenarios_test[2])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_boot_Naive[:,:,1,3],title = "Sampling 2 - method 1, $(scenarios_test[3])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_boot_Heuristic[:,:,1,1],title = "Sampling 2 - method 2, $(scenarios_test[1])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_boot_Heuristic[:,:,1,2],title = "Sampling 2 - method 2, $(scenarios_test[2])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])
pretty_table(gap_boot_Heuristic[:,:,1,3],title = "Sampling 2 - method 2, $(scenarios_test[3])", header = ["$(scenarios[1])", "$(scenarios[2])", "$(scenarios[3])", "$(scenarios[4])", "$(scenarios[5])"])

# Probabilities
prob_gen_Naive = [pro.probability for pro in Stochastic_problem_gen_Naive]
prob_gen_Heuristic = [pro.probability for pro in Stochastic_problem_gen_Heuristic]
prob_boot_Naive = [pro.probability for pro in Stochastic_problem_boot_Naive]
prob_boot_Heuristic = [pro.probability for pro in Stochastic_problem_boot_Heuristic]
p = scatter(prob_gen_Naive[1,1,1,1])
x = [i for i in 1:scenarios[1]]
for i in 1:length(prob_gen_Naive[1,1,1,1])
    annotate!(x[i],prob_gen_Naive[1,1,1,1][i]+0.02, text(prob_gen_Naive[1,1,1,1][i],6))
end
display(p)

prob_gen_Naive[1,1,1,1][1]

########################################################################

# Scenario analysis
plot_folder = "Plots/Data/Scenarios/"
# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end
# Load data - change this to be correct
HPC_folder = "Finlandia_no_cars_light_6014_05_20"
repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
sc = length(scenarios)
n = length(n_unknown)
test_problem = Finlandia_test[5]
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", test_problem, "hazardous"
Deterministic_problem = load_data(problemname1, problemname2, problemname3)
EVP_problem_gen = Array{Any}(undef, repetitions, sc, n)
Stochastic_problem_gen = Array{Any}(undef, repetitions, sc, n)
EVP_problem_boot = Array{Any}(undef, repetitions, sc, n)
Stochastic_problem_boot = Array{Any}(undef, repetitions, sc, n)

for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # Load EVP
            foldername = "EVP_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Problem"
            EVP_problem_gen[i, j, k] = get_deterministic_problem(foldername,
                filename, HPC_folder, problemname1, problemname2, problemname3)
            foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "EVP_Problem"
            EVP_problem_boot[i, j, k] = get_deterministic_problem(foldername,
                filename, HPC_folder, problemname1, problemname2, problemname3)
            # Load Stochastic
            foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem_gen[i, j, k] = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, problemname2, problemname3)
            foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            Stochastic_problem_boot[i, j, k] = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, problemname2, problemname3)
        end
    end
end
# plot for original problem
OG_plot = plot_cargo_OG(Deterministic_problem)
display(OG_plot[1])
display(OG_plot[2])
# Save plots
savefig(OG_plot[1], plot_folder * "Cargo_distribution_Determinstic_1_" * test_problem * ".png")
savefig(OG_plot[2], plot_folder * "Cargo_distribution_Determinstic_2" * test_problem * ".png")
# uniform random sampling
plots_gen = plot_cargo_weights(Stochastic_problem_gen[1, 1, 1], [i for i = 1:scenarios[1]])
xplots = 2
yplots = 2
plot(plots_gen[1], plots_gen[2], plots_gen[3], plots_gen[4], layout=(xplots, yplots))
total_weight_gen = []
weight_diff = []
det_weight = Deterministic_problem.cargo.total_weight
push!(total_weight_gen, Stochastic_problem_gen[1, 1, 1].cargo.items.total_weight)
push!(weight_diff, total_weight_gen[1] .- det_weight)
p = plot([i for i = 1:scenarios[1]], total_weight_gen[1] - ones(scenarios[1]) * det_weight, label="No. scenarios: $(scenarios[1])", xlabel="Scenario", ylabel="Weight (t.)", title="Weight difference")
for i in 2:sc
    push!(total_weight_gen, Stochastic_problem_gen[1, i, 1].cargo.items.total_weight)
    push!(weight_diff, total_weight_gen[i] .- det_weight)
    p = plot!([j for j = 1:scenarios[i]], total_weight_gen[i] - ones(scenarios[i]) * det_weight, label="No. scenarios: $(scenarios[i])")
end
display(p)
savefig(p, plot_folder * "Weight_difference_uniform_random_sampling_nocars_light60.png")
# generate scenarios for other test - should be loaded but is not made yet
gen_problem_nocars_heavy_60 = Array{Any}(undef, sc)
nocars_heavy_60 = Finlandia_test[7]
problemname1, problemname2, problemname3 = "finlandia", nocars_heavy_60, "hazardous"
Det_pro_nocars_heavy_60 = load_data(problemname1, problemname2, problemname3)
det_weight_nocars_heavy_60 = Det_pro_nocars_heavy_60.cargo.total_weight
total_weight_gen_nocars_heavy_60 = []
for i in 1:sc
    gen_problem_nocars_heavy_60[i] = create_stochastic_problem(Det_pro_nocars_heavy_60, scenarios[i], length(Det_pro_nocars_heavy_60.cargo), [])
    push!(total_weight_gen_nocars_heavy_60, gen_problem_nocars_heavy_60[i].cargo.items.total_weight)
end
p = plot([i for i = 1:scenarios[1]], total_weight_gen_nocars_heavy_60[1] - ones(scenarios[1]) * det_weight_nocars_heavy_60, label="No. scenarios: $(scenarios[1])", xlabel="Scenario", ylabel="Weight (t.)", title="Weight difference")
for i in 2:sc
    p = plot!([j for j = 1:scenarios[i]], total_weight_gen_nocars_heavy_60[i] - ones(scenarios[i]) * det_weight_nocars_heavy_60, label="No. scenarios: $(scenarios[i])")
end
display(p)
savefig(p, plot_folder * "Weight_difference_uniform_random_sampling_nocars_heavy60.png")

item_id = 1
one_item_heavy = []
one_item_light = []
for i in 1:scenarios[sc]
    push!(one_item_heavy, gen_problem_nocars_heavy_60[sc].cargo.items[i].items[item_id].weight)
    push!(one_item_light, Stochastic_problem_gen[1, sc, 1].cargo.items[i].items[item_id].weight)
end
p = scatter([1:scenarios[sc]], one_item_heavy, label="Scenario weight", xlabel="Scenario", ylabel="Weight (t.)", title="Weight of cargo id $(item_id)")
p = plot!([1:scenarios[sc]], ones(scenarios[sc]) * Det_pro_nocars_heavy_60.cargo.items[item_id].weight, label="Actual weight")
savefig(p, plot_folder * "Weight_uniform_random_sampling_nocars_heavy60_item_$(item_id).png")
p = scatter([1:scenarios[sc]], one_item_light, label="Scenario weight", xlabel="Scenario", ylabel="Weight (t.)", title="Weight of cargo id $(item_id)")
p = plot!([1:scenarios[sc]], ones(scenarios[sc]) * Deterministic_problem.cargo.items[item_id].weight, label="Actual weight")
savefig(p, plot_folder * "Weight_uniform_random_sampling_nocars_light60_item_$(item_id).png")


# Test scenario reduction

# test speed