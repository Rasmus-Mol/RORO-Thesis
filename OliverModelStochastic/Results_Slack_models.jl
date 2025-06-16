# Script to analyse the results of the Slack models
include("packages_and_files.jl")
# slack results
HPC_folders = [
    "Finlandia_mixed_light_60_SlackDeck1_13_06_18",
    #"Finlandia_mixed_light_100_SlackDeck1_13_06_19",
    "Finlandia_mixed_heavy_60_SlackDeck1_13_06_18",
    #"Finlandia_mixed_heavy_100_SlackDeck1_13_06_18",
    "Finlandia_no_cars_light_60_SlackDeck1_13_06_18",
    #"Finlandia_no_cars_light_100_SlackDeck1_13_06_18",
    "Finlandia_no_cars_heavy_60_SlackDeck1_13_06_18",
    #"Finlandia_no_cars_heavy_100_SlackDeck1_13_06_18",
]
problemname1, problemname3 = "finlandia", "hazardous"
slack_fraction = [[1.1, 1, 1], [1.2, 1, 1], [1.3, 1, 1]]
sl = length(slack_fraction)
repetitions = 5
scenarios = [10]
sc = scenarios[1]
time_limit = 60 * 60

Deterministic_problem = Array{Any}(undef, length(HPC_folders), sl)
Deterministic_Solution = Array{Any}(undef, length(HPC_folders), sl)

EVP_gen_problem = Array{Any}(undef, length(HPC_folders), repetitions, sl)
EVP_gen = Array{Any}(undef, length(HPC_folders), repetitions, sl)
EVP_gen_to_find_mean = Array{Any}(undef, length(HPC_folders), repetitions, sl, sc)
EVP_gen_fitted = Array{Any}(undef, length(HPC_folders), repetitions, sl)
Stochastic_gen = Array{Any}(undef, length(HPC_folders), repetitions, sl)
Stochastic_gen_fitted = Array{Any}(undef, length(HPC_folders), repetitions, sl)

EVP_boot_problem = Array{Any}(undef, length(HPC_folders), repetitions, sl)
EVP_boot = Array{Any}(undef, length(HPC_folders), repetitions, sl)
EVP_boot_to_find_mean = Array{Any}(undef, length(HPC_folders), repetitions, sl, sc)
EVP_boot_fitted = Array{Any}(undef, length(HPC_folders), repetitions, sl)
Stochastic_boot = Array{Any}(undef, length(HPC_folders), repetitions, sl)
Stochastic_boot_fitted = Array{Any}(undef, length(HPC_folders), repetitions, sl)

boot_fitted_slacked = Array{Any}(nothing, length(HPC_folders), repetitions, sl)
gen_fitted_slacked = Array{Any}(nothing, length(HPC_folders), repetitions, sl)
slack_boot = Array{Any}(nothing, length(HPC_folders), repetitions, sl)
slack_gen = Array{Any}(nothing, length(HPC_folders), repetitions, sl)

Stochastic_problem_gen = Array{Any}(undef, length(HPC_folders), repetitions, sl)
Stochastic_problem_boot = Array{Any}(undef, length(HPC_folders), repetitions, sl)

for l in 1:length(HPC_folders)
    test_instance = Finlandia_test[l]
    HPC_folder = HPC_folders[l]
    repetitions_temp, scenarios_temp, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
    #println("Extra info about test: ", note)
    #repetitions = 1 # didn't finish running
    #sc = length(scenarios)
    n = length(n_unknown)
    Deterministic_problem[l] = load_data(problemname1, test_instance, problemname3)
    for j in 1:sl
        Deterministic_Solution[l,j] = get_solution_deterministic("Finlandia_deterministic",
        "Deterministic_Solution_slack_$(j)", HPC_folder)
        for i in 1:repetitions
            # Solutions
            # EVP-Gen
            foldername = "EVP_random_sampling_rep$(i)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)_slack_$(slack_fraction[j])"
            filename = "EVP_Solution"
            #EVP_gen_problem[l, i, j] = get_evp_problem(foldername,
            #    filename, HPC_folder, problemname1, problemname2, problemname3)
            EVP_gen[l,i,j] = get_solution_deterministic(foldername,
            filename,HPC_folder)
            for k in 1:sc
                EVP_gen_to_find_mean[l, i, j, k] = get_solution_deterministic(foldername,
                    "Fitted_Solution_EVP_$(k)", HPC_folder)
            end
            # EVP-boot
            foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)_slack_$(slack_fraction[j])"
            filename = "EVP_Solution"
            #EVP_gen_problem[l, i, j] = get_evp_problem(foldername,
            #    filename, HPC_folder, problemname1, problemname2, problemname3)
            EVP_boot[l,i,j] = get_solution_deterministic(foldername,
            filename,HPC_folder)
            for k in 1:sc
                EVP_boot_to_find_mean[l, i, j, k] = get_solution_deterministic(foldername,
                    "Fitted_Solution_EVP_$(k)", HPC_folder)
            end
            # Stochastic - Gen
            foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)_slack_$(slack_fraction[j])"
            filename = "Stochastic_Problem"
            # problem
            Stochastic_problem_gen[l, i, j] = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, test_instance, problemname3)
            # solution
            filename = "Stochastic_Solution"
            Stochastic_gen[l, i, j] = get_solution_stochastic(foldername,
                filename, HPC_folder)
            Stochastic_gen_fitted[l, i, j] = get_solution_deterministic(foldername,
                "Fitted_Solution", HPC_folder)
            # Load slacked
            if Stochastic_gen_fitted[l, i, j].status != "OPTIMAL"
                gen_fitted_slacked[l, i, j] = get_solution_deterministic(foldername,
                    "Fitted_Solution_slacked", HPC_folder)
                slack_gen[l, i, j] = get_slack(foldername, "Fitted_Solution", HPC_folder)
            end

            # Stochastic Bootstrap
            foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)_slack_$(slack_fraction[j])"
            filename = "Stochastic_Problem"
            # problem
            Stochastic_problem_boot[l, i, j] = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, test_instance, problemname3)
            # solution
            filename = "Stochastic_Solution"
            Stochastic_boot[l, i, j] = get_solution_stochastic(foldername,
                filename, HPC_folder)
            Stochastic_boot_fitted[l, i, j] = get_solution_deterministic(foldername,
                "Fitted_Solution", HPC_folder)
            # Load slacked
            if Stochastic_boot_fitted[l, i, j].status != "OPTIMAL"
                boot_fitted_slacked[l, i, j] = get_solution_deterministic(foldername,
                    "Fitted_Solution_slacked", HPC_folder)
                slack_boot[l, i, j] = get_slack(foldername, "Fitted_Solution", HPC_folder)
            end
        end
    end
end
# how many stowage plans were infeasible
vessel = Deterministic_problem[1,1].vessel
println("Gen:" ,count(!isnothing, gen_fitted_slacked))
println("Boot:" ,count(!isnothing, boot_fitted_slacked))
for i in 1:length(HPC_folders)
    for j in 1:sl
        for k in 1:repetitions
            if sum([c.weight for c in filter(x->x.deck==1 ,Stochastic_gen_fitted[i, k, j].cargo)]) > vessel.decks[1].weight_limit*slack_fraction[j][1]
                println("Stochastic_gen_fitted $(i), $(j), $(k) has too much cargo on deck 1")
            end
            if sum([c.weight for c in filter(x-> x.deck == 1,Stochastic_boot_fitted[i, k, j].cargo)]) > vessel.decks[1].weight_limit*slack_fraction[j][1]
                println("Stochastic_boot_fitted $(i), $(j), $(k) has too much cargo on deck 1")
            end
        end
    end
end
vessel.decks[1].weight_limit*slack_fraction[1][1]
sum([c.weight for c in filter(x->x.deck==1 ,Stochastic_gen_fitted[i, k, j].cargo)])
#
sum([c.weight for c in filter(x-> x.deck == 1,Stochastic_boot_fitted[4,1,1].cargo)])


cargo_loaded_gen = Array{Any}(undef, length(HPC_folders), repetitions, sl)
cargo_loaded_boot = Array{Any}(undef, length(HPC_folders), repetitions, sl)
objective_val = Array{Any}(undef, length(HPC_folders), repetitions, sl)
objective_val_boot = Array{Any}(undef, length(HPC_folders), repetitions, sl)
gap_gen = Array{Any}(undef, length(HPC_folders), repetitions, sl)
gap_boot = Array{Any}(undef, length(HPC_folders), repetitions, sl)
for i in 1:length(HPC_folders)
    # save cargo Loaded
    for j in 1:sl
        for k in 1:repetitions
            cargo_loaded_gen[i, k, j] = Stochastic_gen[i, k, j].n_cargo_loaded
            cargo_loaded_boot[i, k, j] = Stochastic_boot[i, k, j].n_cargo_loaded
            objective_val[i, k, j] = Stochastic_gen[i, k, j].objective
            objective_val_boot[i, k, j] = Stochastic_boot[i, k, j].objective
            gap_gen[i, k, j] = round(Stochastic_gen[i, k, j].gap *100, digits = 2)
            gap_boot[i, k, j] = round(Stochastic_boot[i, k, j].gap*100,digits=2)
        end
    end

end
avg_cargo_loaded = dropdims(sum(cargo_loaded_gen, dims=2) ./ repetitions; dims =2)
avg_cargo_loaded_boot = dropdims(sum(cargo_loaded_boot, dims=2) ./ repetitions; dims = 2)

gap_gen[1,:,1]
objective_val_boot[1,:,1]


# EVP
temp = mean(getfield.(EVP_boot_to_find_mean[1,1,1,:],:objective))
Stochastic_boot[1,1,1].objective
println(getfield.(EVP_boot_to_find_mean[1,1,1,:],:objective))
println(Stochastic_boot[1,1,1].objective)

Deterministic_Solution[1,1].objective
println(Stochastic_boot_fitted[4,1,1].objective)
println(Stochastic_gen_fitted[4,1,1].objective)


