# Main for stability test for stochastic solutions
include("packages_and_files.jl")

HPC_folder = "Finlandia_stability_tests_stochastic"
problemname1, problemname3 = "finlandia", "hazardous"
no_test = 10
repetitions = 1#10  # number of repetitions
sc = 5 # number of scenarios
scenarios = [10,20,30,40,50]

# Test problems and solutions
# A problem for each instance for each repetition run
Sto_pro_gen = Array{Any}(undef, 8, sc, repetitions)
slack_sol_gen = Array{Any}(nothing,8, sc, repetitions, no_test)
sol_gen = Array{Any}(undef, 8, sc, repetitions, no_test)

Sto_pro_boot = Array{Any}(undef, 8, sc, repetitions)
slack_sol_boot = Array{Any}(nothing,8, sc, repetitions, no_test)
sol_boot = Array{Any}(undef, 8, sc, repetitions, no_test)

for i in 1:8
    test_problem_name = Finlandia_test[i]
    for l in 1:sc
        for j in 1:repetitions
            # Load the problem
            foldername = "Stochastic_Stability_randomscenarios_uniformsampling_instance_$(i)_scenarios_$(scenarios[l])_repetition_$(j)"
            filename = "Stochastic_Problem"
            Sto_pro_gen[i, l, j] = get_stochastic_problem(foldername, filename, HPC_folder,
                problemname1, test_problem_name, problemname3)
            # Load the solution
            path = joinpath("Results", HPC_folder, foldername)
            files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
            for k in 1:no_test
                if "Solution_$(k)_info.json" in files # Could be solved
                    sol_gen[i, l, j, k] = get_solution_deterministic(foldername, "Solution_$(k)", HPC_folder)
                else # was infeasible and therefore slacked
                    sol_gen[i, l, j, k] = get_solution_deterministic(foldername, "Solution_$(k)_slacked", HPC_folder)
                    # Should work when next test has been run 
                    slack_sol_gen[i, l, j, k] = get_slack(foldername, "Fitted_Solution_slacked_$(k)", HPC_folder)
                    #slack_sol_gen_idx[i,j] = 1
                end
            end
            # Load the problem
            foldername = "Stochastic_Stability_randomscenarios_Bootstrap_instance_$(i)_scenarios_$(scenarios[l])_repetition_$(j)"
            filename = "Stochastic_Problem"
            Sto_pro_boot[i, l, j] = get_stochastic_problem(foldername, filename, HPC_folder,
                problemname1, test_problem_name, problemname3)
            # Load the solution
            path = joinpath("Results", HPC_folder, foldername)
            files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
            for k in 1:no_test
                if "Solution_$(k)_info.json" in files # Could be solved
                    sol_boot[i, l, j, k] = get_solution_deterministic(foldername, "Solution_$(k)", HPC_folder)
                else # was infeasible and therefore slacked
                    sol_boot[i, l, j, k] = get_solution_deterministic(foldername, "Solution_$(k)_slacked", HPC_folder)
                    # Should work when next test has been run 
                    slack_sol_boot[i, l, j, k] = get_slack(foldername, "Fitted_Solution_slacked_$(k)", HPC_folder)
                    #slack_sol_gen_idx[i,j] = 1
                end
            end
        end
    end
end

inf_count_gen = zeros(8, sc, repetitions)
inf_count_boot = zeros(8, sc, repetitions)
# count how many of the test were feasible
for i in 1:8 # number of instances
    for l in 1:sc
        for j in 1:repetitions # number of repetitions
            for k in 1:no_test
                if !isnothing(slack_sol_gen[i, l, j, k]) # if there is a slack solution
                    inf_count_gen[i, l, j] += 1 # count as infeasible
                    println("instance: $(i), Repetition: $(j), Test: $(k) was infeasible")
                    if sum(slack_sol_gen[i, l, j, k][1]) < 1 # if there is a slack
                        println("Not deck 1 weight limit that is the issue: ", slack_sol_gen[i, l, j, k])
                    end
                end
                if !isnothing(slack_sol_boot[i, l, j, k]) # if there is a slack solution
                    inf_count_boot[i, l, j] += 1 # count as infeasible
                    println("instance: $(i), Repetition: $(j), Test: $(k) was infeasible")
                    if sum(slack_sol_boot[i, l, j, k][1]) < 1 # if there is a slack
                        println("Not deck 1 weight limit that is the issue: ", slack_sol_boot[i, l, j, k])
                    end
                end
            end
        end
    end
end
inf_count_gen[:,1,1]
inf_count_gen
inf_count_boot 

slack_sol_gen[8,1,1,1]
sol_gen[1,1,1,1].ballast_weight
sol_gen[1,1,1,4].ballast_weight

sum(sol_gen[8,1,1,1].cs)
# OOPS sol_gen have the deterministic cargo weight I think
for i in 1:no_test
    #println(sum([c.weight for c in filter(x->x.deck == 1, sol_gen[8,1,1,i].cargo)]))
    if !isnothing(slack_sol_gen[8,1,1,i])
        println("Slack: ", sum(slack_sol_gen[8,1,1,i][1]))
    end
end
temp = sol_gen[8,1,1,1].cs
sol_gen[8,1,1,1].cargo_weight
sum(temp)
findall(x-> sum(temp[x,:]) < 0.5 ,1:size(temp, 1))
Sto_pro_gen[8,1,1].cargo.


Sto_pro_gen[8,1,1].cargo.items.total_weight

# Change in ballast water
ballast_water_gen = zeros(8,no_test)
ballast_water_boot = zeros(8,no_test)
avg_ballast_water_gen = zeros(8)
avg_ballast_water_boot = zeros(8)
for i in 1:8
    for j in 1:no_test
        if isnothing(slack_sol_gen[i, 1, 1, j])# was feasible
            ballast_water_gen[i,j] = sol_gen[i,1,1,j].ballast_weight
        end
        if isnothing(slack_sol_boot[i, 1, 1, j])# was feasible
            ballast_water_boot[i,j] = sol_boot[i,1,1,j].ballast_weight
        end
    end
    println("####################")
    #println("Deterministic ballast water: ", Det_sol[i].ballast_weight)
    println("Ballast water gen: ", ballast_water_gen[i,:])
    println("Ballast water boot: ", ballast_water_boot[i,:])
    avg_ballast_water_boot[i] = mean(filter(x-> x>0,ballast_water_boot[i,:]))
    avg_ballast_water_gen[i] = mean(filter(x-> x>0,ballast_water_gen[i,:]))
end



# Stochastic solution results

HPC_folders = [
    "Finlandia_mixed_light_60_28_05_13",
    "Finlandia_mixed_light_100_28_05_13",
    "Finlandia_mixed_heavy_60_28_05_14",
    "Finlandia_mixed_heavy_100_28_05_15",
    "Finlandia_no_cars_light_60_28_05_16",
    "Finlandia_no_cars_light_100_28_05_16",
    "Finlandia_no_cars_heavy_60_28_05_19",
    "Finlandia_no_cars_heavy_100_28_05_19",
]

# HPC_folder
test_instance = Finlandia_test[1]
HPC_folder = HPC_folders[1]#"Finlandia_"*test_instance*"_15_05_09"
#=
plot_folder = "Plots/Results/Finlandia_"*test_instance*"/"
# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end=#
# Load data from HPC
repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
println("Extra info about test: ", note)
#repetitions = 1 # didn't finish running
sc = length(scenarios)
n = length(n_unknown)
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", test_instance, "hazardous"
Deterministic_problem = Array{Any}(undef, length(HPC_folders))
#load_data(problemname1, problemname2, problemname3)
Deterministic_Solution = Array{Any}(undef, length(HPC_folders))
#get_solution_deterministic("Finlandia_deterministic",
#    "Deterministic_Solution", HPC_folder)
# soluton arrays
#EVP_gen = Array{Any}(undef, repetitions, sc,n)
#EVP_gen_fitted = Array{Any}(undef, repetitions, sc,n)
Stochastic_gen = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
Stochastic_gen_fitted = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
#EVP_boot = Array{Any}(undef, repetitions, sc,n)
#EVP_boot_fitted = Array{Any}(undef, repetitions, sc,n)
Stochastic_boot = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
Stochastic_boot_fitted = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))

boot_fitted_slacked = Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))
gen_fitted_slacked = Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))
slack_boot= Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))
slack_gen = Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))

Stochastic_problem_gen = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
Stochastic_problem_boot = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))

for l in 1:length(HPC_folders)
        test_instance = Finlandia_test[l]
        HPC_folder = HPC_folders[l]
        repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
        println("Extra info about test: ", note)
        #repetitions = 1 # didn't finish running
        #sc = length(scenarios)
        n = length(n_unknown)
        Deterministic_Solution[l] = get_solution_deterministic("Finlandia_deterministic",
            "Deterministic_Solution", HPC_folder)
        Deterministic_problem[l] = load_data(problemname1, test_instance, problemname3)
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
                filename = "Stochastic_Problem"
                # problem
                Stochastic_problem_gen[i,j,k,l] = get_stochastic_problem(foldername,
                    filename,HPC_folders[l],problemname1,problemname2,problemname3)
                # solution
                filename = "Stochastic_Solution"
                Stochastic_gen[i, j, k,l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
                Stochastic_gen_fitted[i, j, k,l] = get_solution_deterministic(foldername,
                    "Fitted_Solution", HPC_folder)
                # Load slacked
                if Stochastic_gen_fitted[i, j, k,l].status != "OPTIMAL"
                    gen_fitted_slacked[i,j,k,l] = get_solution_deterministic(foldername,
                    "Fitted_Solution_slacked", HPC_folder)
                    slack_gen[i,j,k,l] = get_slack(foldername, "Fitted_Solution",HPC_folder)
        
                end
                # Stochastic Bootstrap
                foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
                filename = "Stochastic_Problem"
                # problem
                Stochastic_problem_boot[i,j,k,l] = get_stochastic_problem(foldername,
                    filename,HPC_folders[l],problemname1,problemname2,problemname3)
                # solution
                filename = "Stochastic_Solution"
                Stochastic_boot[i, j, k,l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
                Stochastic_boot_fitted[i, j, k,l] = get_solution_deterministic(foldername,
                    "Fitted_Solution", HPC_folder)
                # Load slacked
                if Stochastic_boot_fitted[i, j, k,l].status != "OPTIMAL"
                    boot_fitted_slacked[i,j,k,l] = get_solution_deterministic(foldername,
                    "Fitted_Solution_slacked", HPC_folder)
                    slack_boot[i,j,k,l] = get_slack(foldername, "Fitted_Solution",HPC_folder)
                end
            end
        end
    end
end