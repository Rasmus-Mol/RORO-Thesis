# Main script to test stability of stochastic solutions
include("packages_and_files.jl")
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
slack_boot = Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))
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
                Stochastic_problem_gen[i, j, k, l] = get_stochastic_problem(foldername,
                    filename, HPC_folders[l], problemname1, problemname2, problemname3)
                # solution
                filename = "Stochastic_Solution"
                Stochastic_gen[i, j, k, l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
                Stochastic_gen_fitted[i, j, k, l] = get_solution_deterministic(foldername,
                    "Fitted_Solution", HPC_folder)
                # Load slacked
                if Stochastic_gen_fitted[i, j, k, l].status != "OPTIMAL"
                    gen_fitted_slacked[i, j, k, l] = get_solution_deterministic(foldername,
                        "Fitted_Solution_slacked", HPC_folder)
                    slack_gen[i, j, k, l] = get_slack(foldername, "Fitted_Solution", HPC_folder)

                end
                # Stochastic Bootstrap
                foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
                filename = "Stochastic_Problem"
                # problem
                Stochastic_problem_boot[i, j, k, l] = get_stochastic_problem(foldername,
                    filename, HPC_folders[l], problemname1, problemname2, problemname3)
                # solution
                filename = "Stochastic_Solution"
                Stochastic_boot[i, j, k, l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
                Stochastic_boot_fitted[i, j, k, l] = get_solution_deterministic(foldername,
                    "Fitted_Solution", HPC_folder)
                # Load slacked
                if Stochastic_boot_fitted[i, j, k, l].status != "OPTIMAL"
                    boot_fitted_slacked[i, j, k, l] = get_solution_deterministic(foldername,
                        "Fitted_Solution_slacked", HPC_folder)
                    slack_boot[i, j, k, l] = get_slack(foldername, "Fitted_Solution", HPC_folder)
                end
            end
        end
    end
end

# Generate random tests and see if the solutions are stable
HPC_folder_save = "Finlandia_stability_tests_stochastic/"
if !isdir("Results/" * HPC_folder_save)
    mkdir("Results/" * HPC_folder_save)
end

no_test = 10

for i in 1:length(HPC_folders)
    test_problem_name = Finlandia_test[i]
    HPC_folder_load = HPC_folders[i]
    # Load solution
    det_sol = Deterministic_Solution[i]
    det_pro = Deterministic_problem[i]
    slots = det_pro.slots
    vessel = det_pro.vessel
    for l in 1:sc
        #for k in 1:repetitions
            sto_sol = Stochastic_boot[1, l, 1, i]
            # 
            sto_pro = create_stochastic_problem(det_pro, no_test, length(det_pro.cargo), [])
            # test
            #HPC_folder_save_temp = HPC_folder_save*test_problem_name*"/" # testing all problems
            foldername = "Stochastic_Stability_randomscenarios_uniformsampling_instance_$(i)_scenarios_$(scenarios[l])_repetition_$(1)"#repetition_$(k)"
            # save scenario
            filename = "Stochastic_Problem"
            write_problem_stochastic(sto_pro, foldername, filename, HPC_folder_save)
            for j in 1:no_test
                new_c = sto_pro.cargo.items[j]
                f, strs, mo = feasibility_check_stochastic(sto_sol, det_pro, new_c)
                if f == true # found solution
                    sol = get_solution_second_stage_stochastic(det_pro, mo, sto_sol)
                    # save result
                    filename = "Solution_$(j)"
                    write_solution(sol, foldername, filename, HPC_folder_save)
                else # Problem was infeasible
                    fitted_sol_slacked = get_solution_second_stage_stochastic(det_pro, mo, sto_sol)
                    #get_solution_second_stage_deterministic(det_pro, mo, det_sol)
                    write_solution(fitted_sol_slacked, foldername, "Solution_$(j)_slacked", HPC_folder_save)
                    write_slack(HPC_folder_save, foldername, "Fitted_Solution_slacked_$(j)", mo)
                    #infeasible_test1_gen[j,i] = 1 # infeasible
                end
            end
            # Other sceario method
            sto_pro = create_stochastic_problem(det_pro, no_test, length(det_pro.cargo), [],Bootstrap_bookedweight_quantile) 
            # test
            #HPC_folder_save_temp = HPC_folder_save*test_problem_name # testing all problems
            foldername = "Stochastic_Stability_randomscenarios_Bootstrap_instance_$(i)_scenarios_$(scenarios[l])_repetition_$(1)"#$(k)"
            # save scenario
            filename = "Stochastic_Problem"
            write_problem_stochastic(sto_pro, foldername, filename, HPC_folder_save)
            for j in 1:no_test
                new_c = sto_pro.cargo.items[j]
                f, strs, mo = feasibility_check_stochastic(sto_sol, det_pro, new_c)
                if f == true # found solution
                    sol = get_solution_second_stage_stochastic(det_pro, mo, sto_sol)
                    # save result
                    filename = "Solution_$(j)"
                    write_solution(sol, foldername, filename, HPC_folder_save)
                else # Problem was infeasible
                    fitted_sol_slacked = get_solution_second_stage_stochastic(det_pro, mo, sto_sol)
                    #get_solution_second_stage_deterministic(det_pro, mo, det_sol)
                    write_solution(fitted_sol_slacked, foldername, "Solution_$(j)_slacked", HPC_folder_save)
                    write_slack(HPC_folder_save, foldername, "Fitted_Solution_slacked_$(j)", mo)
                    #infeasible_test1_gen[j,i] = 1 # infeasible
                end
            end   
        #end
    end
end