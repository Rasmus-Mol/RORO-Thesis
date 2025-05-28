include("packages_and_files.jl")

problemname1, problemname3 = "finlandia", "hazardous"
HPC_folders = [
    "Finlandia_mixed_light_60_15_05_13",
    "Finlandia_mixed_light_100_15_05_15",
    "Finlandia_mixed_heavy_60_15_05_17",
    "Finlandia_mixed_heavy_100_15_05_17",
    "Finlandia_no_cars_light_60_14_05_20",
    "Finlandia_no_cars_light_100_15_05_10",
    "Finlandia_no_cars_heavy_60_15_05_10",
    "Finlandia_no_cars_heavy_100_15_05_09",
]
no_test1 = 10
no_test2 = 10

# Load problems and results
no_folders = length(HPC_folders)
foldername = "Determinitic_Stability_randomscenarios_uniformsampling"
test1_problems_gen = Array{Any}(undef, no_folders)
test1_problems_boot = Array{Any}(undef, no_folders)
test1_solutions_gen = Array{Any}(undef, no_folders,no_test1)
test1_solutions_boot = Array{Any}(undef, no_folders,no_test1)
#slack_sol_gen = Array{Any}(undef, no_folders,no_test1)
#slack_sol_boot = Array{Any}(undef, no_folders,no_test1)
slack_sol_gen = zeros(no_folders,test1)
slack_sol_boot = zeros(no_folders,test1)

initial_random_plan_trucksfirst = Array{Any}(undef, no_folders,no_test2)
initial_random_plan_carsfirst = Array{Any}(undef, no_folders,no_test2)
initial_random_plan_secufirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_trucksfirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_carsfirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_secufirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_trucksfirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_carsfirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_secufirst = Array{Any}(undef, no_folders,no_test2)

for i in 1:1#no_folders
    test_problem_name = Finlandia_test[i]
    HPC_folder_load = HPC_folders[i]
    # Test 1 - Uniform sampling
    foldername = "Determinitic_Stability_randomscenarios_uniformsampling"
    filename = "Stochastic_Problem"
    test1_problems_gen[i] = get_stochastic_problem(foldername,
                filename, HPC_folder_load, problemname1, test_problem_name, problemname3)
    path = joinpath("Results", HPC_folder_load,foldername)
    files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
    for j in 1:no_test1
        if "Solution_$(i)_info.json" in files # Could be solved
            test1_solutions_gen[i,j] = get_solution_deterministic(foldername,"Solution_$(i)",HPC_folder_load)
        else # was infeasible and therefore slacked
            test1_solutions_gen[i,j] = get_solution_deterministic(foldername,"Solution_$(i)_slacked",HPC_folder_load)
            # Should work when next test has been run 
            #slack_sol_gen[i,j]  = get_slack(foldername,"Fitted_Solution_slacked_$(i)",HPC_folder_load)
            slack_sol_gen[i,j] = 1
        end
    end
    foldername = "Determinitic_Stability_randomscenarios_bootstrapsampling"
    filename = "Stochastic_Problem"
    test1_problems_boot[i] = get_stochastic_problem(foldername,
                filename, HPC_folder_load, problemname1, test_problem_name, problemname3)
    path = joinpath("Results", HPC_folder_load,foldername)
    files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
    for j in 1:no_test1
        if "Solution_$(i)_info.json" in files # Could be solved
            test1_solutions_boot[i,j] = get_solution_deterministic(foldername,"Solution_$(i)",HPC_folder_load)
        else # was infeasible and therefore slacked
            test1_solutions_boot[i,j] = get_solution_deterministic(foldername,"Solution_$(i)_slacked",HPC_folder_load)
            # Should work when next test has been run 
            #slack_sol_boot[i,j] = get_slack(foldername,"Fitted_Solution_slacked_$(i)",HPC_folder_load)
            slack_sol_boot[i,j] = 1
        end
    end
    # Test 2 - random plan
    for j in 1:no_test2
        # Trucks first
        foldername = "Random_Stowage_Plan_Trucks_first"
        filename = "Random_plan" # should be fixed next time
        # filename = "Random_plan_$(j)" 
        initial_random_plan_trucksfirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_trucksfirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(i)",HPC_folder_load)
        sol_shifts_trucksfirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(i)",HPC_folder_load)
        # Cars first
        foldername = "Random_Stowage_Plan_Cars_first"
        filename = "Random_plan" # should be fixed next time
        # filename = "Random_plan_$(j)" 
        initial_random_plan_carsfirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_carsfirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(i)",HPC_folder_load)
        sol_shifts_carsfirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(i)",HPC_folder_load)
        # Secu first
        foldername = "Random_Stowage_Plan_Secu_first"
        filename = "Random_plan" # should be fixed next time
        # filename = "Random_plan_$(j)" 
        initial_random_plan_secufirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_secufirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(i)",HPC_folder_load)
        sol_shifts_secufirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(i)",HPC_folder_load) 
    end

end
# Empty ship
HPC_folder_load = HPC_folders[1]
sol_empty_ship = get_solution_deterministic("Finlandia_deterministic","Empty_ship", HPC_folder_load)