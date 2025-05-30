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
test1_problems_gen = Array{Any}(undef, no_folders)
test1_problems_boot = Array{Any}(undef, no_folders)
test1_solutions_gen = Array{Any}(undef, no_folders,no_test1)
test1_solutions_boot = Array{Any}(undef, no_folders,no_test1)
slack_sol_gen = Array{Any}(undef, no_folders,no_test1)
slack_sol_boot = Array{Any}(undef, no_folders,no_test1)
slack_sol_gen_idx = zeros(no_folders,no_test1)
slack_sol_boot_idx = zeros(no_folders,no_test1)

initial_random_plan_trucksfirst = Array{Any}(undef, no_folders,no_test2)
initial_random_plan_carsfirst = Array{Any}(undef, no_folders,no_test2)
initial_random_plan_secufirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_trucksfirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_carsfirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_secufirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_trucksfirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_carsfirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_secufirst = Array{Any}(undef, no_folders,no_test2)

Det_sol = Array{Any}(undef, no_folders)
Det_pro = Array{Any}(undef, no_folders)

for i in 1:no_folders
    test_problem_name = Finlandia_test[i]
    HPC_folder_load = HPC_folders[i]
    # load deterministic solution and problem
    Det_sol[i] =  get_solution_deterministic("Finlandia_deterministic",
    "Deterministic_Solution", HPC_folder_load)
    Det_pro[i] = load_data(problemname1, test_problem_name, problemname3)
    # Test 1 - Uniform sampling
    foldername = "Determinitic_Stability_randomscenarios_uniformsampling"
    filename = "Stochastic_Problem"
    test1_problems_gen[i] = get_stochastic_problem(foldername,
                filename, HPC_folder_load, problemname1, test_problem_name, problemname3)
    path = joinpath("Results", HPC_folder_load,foldername)
    files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
    for j in 1:no_test1
        if "Solution_$(j)_info.json" in files # Could be solved
            test1_solutions_gen[i,j] = get_solution_deterministic(foldername,"Solution_$(j)",HPC_folder_load)
        else # was infeasible and therefore slacked
            test1_solutions_gen[i,j] = get_solution_deterministic(foldername,"Solution_$(j)_slacked",HPC_folder_load)
            # Should work when next test has been run 
            slack_sol_gen[i,j]  = get_slack(foldername,"Fitted_Solution_slacked_$(j)",HPC_folder_load)
            slack_sol_gen_idx[i,j] = 1
        end
    end
    foldername = "Determinitic_Stability_randomscenarios_bootstrapsampling"
    filename = "Stochastic_Problem"
    test1_problems_boot[i] = get_stochastic_problem(foldername,
                filename, HPC_folder_load, problemname1, test_problem_name, problemname3)
    path = joinpath("Results", HPC_folder_load,foldername)
    files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
    for j in 1:no_test1
        if "Solution_$(j)_info.json" in files # Could be solved
            test1_solutions_boot[i,j] = get_solution_deterministic(foldername,"Solution_$(j)",HPC_folder_load)
        else # was infeasible and therefore slacked
            test1_solutions_boot[i,j] = get_solution_deterministic(foldername,"Solution_$(j)_slacked",HPC_folder_load)
            # Should work when next test has been run 
            slack_sol_boot[i,j] = get_slack(foldername,"Fitted_Solution_slacked_$(i)",HPC_folder_load)
            slack_sol_boot_idx[i,j] = 1
        end
    end
    # Test 2 - random plan
    for j in 1:no_test2
        # Trucks first
        foldername = "Random_Stowage_Plan_Trucks_first"
        #filename = "Random_plan" # should be fixed next time
        filename = "Random_plan_$(j)" 
        initial_random_plan_trucksfirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_trucksfirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(j)",HPC_folder_load)
        sol_shifts_trucksfirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(j)",HPC_folder_load)
        # Cars first
        foldername = "Random_Stowage_Plan_Cars_first"
        #filename = "Random_plan" # should be fixed next time
        filename = "Random_plan_$(j)" 
        initial_random_plan_carsfirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_carsfirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(j)",HPC_folder_load)
        sol_shifts_carsfirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(j)",HPC_folder_load)
        # Secu first
        foldername = "Random_Stowage_Plan_Secu_first"
        #filename = "Random_plan" # should be fixed next time
        filename = "Random_plan_$(j)" 
        initial_random_plan_secufirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_secufirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(j)",HPC_folder_load)
        sol_shifts_secufirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(j)",HPC_folder_load) 
    end

end
# Empty ship
HPC_folder_load = HPC_folders[1]
sol_empty_ship = get_solution_deterministic("Finlandia_deterministic","Empty_ship", HPC_folder_load)
# Empty ship and no ballast water
empty_ship_no_ballast = empty_ship_model(Det_pro[1])
optimize!(empty_ship_no_ballast)
# check if model is infeasible
if termination_status(empty_ship_no_ballast) in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
    println("Empty ship model is infeasible")
else
    println("Empty ship model is feasible")
end
#######################
# test 1
HPC_folders[6]
idx_gen = findall(x -> x==1, slack_sol_gen_idx)
idx_boot = findall(x -> x==1, slack_sol_boot_idx)
println("Test that were infeasible - gen: ",idx_gen)
println("Test that were infeasible - boot: ",idx_boot)
println("number of test that were unfeasible: ", Int(sum(slack_sol_gen_idx)),"/", no_folders*no_test1)
println("number of test that were unfeasible boot: ", Int(sum(slack_sol_boot_idx)),"/", no_folders*no_test1)
# Problems that were still feasible
for i in 1:no_folders
    println("Deterministic total weight: ",Det_pro[i].cargo.total_weight)
    for j in 1:no_test1
        if slack_sol_gen_idx[i,j] == 0 # was feasible
            println("scenario total weight: ",test1_problems_gen[i].cargo.items[j].total_weight)
        end
    end
end
# Change in ballast water
ballast_water_gen = zeros(no_folders,no_test1)
ballast_water_boot = zeros(no_folders,no_test1)
avg_ballast_water_gen = zeros(no_folders)
avg_ballast_water_boot = zeros(no_folders)
for i in 1:no_folders
    for j in 1:no_test1
        if slack_sol_gen_idx[i,j] == 0 # was feasible
            ballast_water_gen[i,j] = test1_solutions_gen[i,j].ballast_weight
        end
        if slack_sol_boot_idx[i,j] == 0 # was feasible
            ballast_water_boot[i,j] = test1_solutions_boot[i,j].ballast_weight
        end
    end
    println("####################")
    println("Deterministic ballast water: ", Det_sol[i].ballast_weight)
    println("Ballast water gen: ", ballast_water_gen[i,:])
    println("Ballast water boot: ", ballast_water_boot[i,:])
    avg_ballast_water_boot[i] = mean(filter(x-> x>0,ballast_water_boot[i,:]))
    avg_ballast_water_gen[i] = mean(filter(x-> x>0,ballast_water_gen[i,:]))
end
for i in 1:no_folders
    #println(round(Det_sol[i].ballast_weight,digits = 2))
    #println(round(Det_pro[i].cargo.total_weight,digits = 2))
    #println(sum(Det_pro[i].cargo.items[j].))
    println(round(mean(test1_problems_boot[i].cargo.items.total_weight),digits = 2))
    #println(no_test1-sum(slack_sol_gen_idx[i,:]))
    #println(round(avg_ballast_water_gen[i],digits=2))
end


test1_problems_boot[1].cargo.items.total_weight
test1_problems_gen[1].cargo.items.total_weight

# infeasible because: Needs new results before doing
inf_problems_gen = [(sum(slack_sol_gen_idx[i,:])>0) ? i : 0 for i in 1:no_folders]
inf_problem_boot = [(sum(slack_sol_boot_idx[i,:])>0) ? i : 0 for i in 1:no_folders]
for i in 1:no_folders
    println("Deterministic total weight: ",Det_pro[i].cargo.total_weight)
    for j in 1:no_test1
        if slack_sol_gen_idx[i,j] == 1 # was infeasible
            println("scenario total weight: ",test1_solutions_gen[i,j].cargo.total_weight)
        end
    end
end
# Slack order:
# slack_deck, slack_Vmax, slack_Vmin, slack_Tmin, slack_Tmax,
# slack_Lmin, slack_Lmax, slack_shearMin, slack_shearMax, slack_bendingMax,
# slack_ballast_tanks





#################################

# test 2 - Stability of random plan
shifts_carsfirst = zeros(no_folders,no_test2)
shifts_trucksfirst = zeros(no_folders,no_test2)
shifts_secufirst = zeros(no_folders,no_test2)
for i in 1:no_folders
    for j in 1:no_test2
        shifts_trucksfirst[i,j] = count((!=).(initial_random_plan_trucksfirst[i,j].cs,
                                sol_no_shifts_trucksfirst[i,j].cs))
        shifts_carsfirst[i,j] = count((!=).(initial_random_plan_carsfirst[i,j].cs,
                                sol_no_shifts_carsfirst[i,j].cs))
        shifts_secufirst[i,j] = count((!=).(initial_random_plan_secufirst[i,j].cs,
                                sol_no_shifts_secufirst[i,j].cs))
    end 
end
# Right now only last column is correct
println(shifts_trucksfirst[:,end])
println(shifts_carsfirst[:,end])
println(shifts_secufirst[:,end])



# Vi har random plan cs
# no shifts solution
# shift allowed solution

