# Script to run on HPC
#=
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
using .StowagePlannerStochastic
using JuMP
=#

include("packages_and_files.jl")

#parse_index = parse(Int, ARGS[1]) # Jobindex input
# Choose instance:

# test 1 - load solution and change weight and test feasibility.
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
#HPC_folders = ["Finlandia_mixed_light_60_15_05_13"]
tests = length(HPC_folders)
no_test = 10
no_test2 = 10

infeasible_test1_gen = zeros(tests,no_test)
infeasible_test2_boot = zeros(tests,no_test)

for j in 1:length(HPC_folders)
    println("Iteration $(j) of $(length(HPC_folders))")
    test_problem_name = Finlandia_test[j]
    HPC_folder_load = HPC_folders[j]
    # Load solution
    det_sol = get_solution_deterministic("Finlandia_deterministic",
        "Deterministic_Solution", HPC_folder_load)
    det_pro = load_data(problemname1, test_problem_name, problemname3)
    slots = det_pro.slots
    vessel = det_pro.vessel
    #=
    sto_pro = create_stochastic_problem(det_pro, no_test, length(det_pro.cargo), [])
    # test
    HPC_folder_save = HPC_folder_load # save in same folder
    foldername = "Determinitic_Stability_randomscenarios_uniformsampling"
    # save scenario
    filename = "Stochastic_Problem"
    write_problem_stochastic(sto_pro, foldername, filename, HPC_folder_save)
    # Test 1
    for i in 1:no_test
        new_c = sto_pro.cargo.items[i]
        f, strs, mo = feasibility_check(det_sol, det_pro, new_c)
        if f == true # found solution
            sol = get_solution_second_stage_deterministic(det_pro, mo, det_sol)
            # save result
            filename = "Solution_$(i)"
            write_solution(sol, foldername, filename, HPC_folder_save)
        else # Problem was infeasible
            fitted_sol_slacked = get_solution_second_stage_deterministic(det_pro, mo, det_sol)
            write_solution(fitted_sol_slacked,foldername,"Solution_$(i)_slacked",HPC_folder_save)
            write_slack(HPC_folder_save, foldername, "Fitted_Solution_slacked_$(i)", mo)
            #infeasible_test1_gen[j,i] = 1 # infeasible
        end
    end
    =#
    # Now using other method for generating scenarios
    sto_pro = create_stochastic_problem(det_pro, no_test, length(det_pro.cargo), [],Bootstrap_bookedweight_quantile) 
    # test
    HPC_folder_save = HPC_folder_load # save in same folder
    foldername = "Determinitic_Stability_randomscenarios_bootstrapsampling"
    # save scenario
    filename = "Stochastic_Problem"
    write_problem_stochastic(sto_pro, foldername, filename, HPC_folder_save)
    # Test 1
    for i in 1:no_test
        new_c = sto_pro.cargo.items[i]
        f, strs, mo = feasibility_check(det_sol, det_pro, new_c)
        if f == true # found solution
            sol = get_solution_second_stage_deterministic(det_pro, mo, det_sol)
            # save result
            filename = "Solution_$(i)"
            write_solution(sol, foldername, filename, HPC_folder_save)
        else
            fitted_sol_slacked = get_solution_second_stage_deterministic(det_pro, mo, det_sol)
            write_solution(fitted_sol_slacked,foldername,"Solution_$(i)_slacked",HPC_folder_save)
            write_slack(HPC_folder_save, foldername, "Fitted_Solution_slacked_$(i)", mo)
            infeasible_test2_boot[j,i] = 1 # infeasible
        end
    end
end
    #####################
    #=
    # Test 2 - Random stowage plan
    cargoc = det_pro.cargo
    ntypes = [length(findall(x->x.cargo_type_id == i, cargoc)) for i in 1:4]
    # Random stowage plan, placed in cargoc order ie cargoc order: 
    # truck, car, Machine, Secu
    foldername = "Random_Stowage_Plan_Trucks_first"
    for i in 1:no_test
        cs_random, not_stowed = random_stowage_plan(cargoc, slots)
        random_plan_not_solved = Solution(
            gap = Inf,
            status = 0,
            objective = Inf,
            time = 0,
            cargo_weight = 0.0,
            total_weight = 0.0,
            ballast_weight = 0.0, 
            area_utilization = 0.0,
            cargo = CargoPlacement[],
            cs = cs_random,
            n_cargo_total = length(cargoc),
            n_cargo_loaded = length(cargoc)-length(not_stowed),
            shear_force = Float64[],
            bending_moment = Float64[],
            ballast_volume = Float64[],
            lcg = 0.0,
            tcg = 0.0, 
            vcg = 0.0,
            n_variables = 0,
            n_binary_variables = 0,
            n_constraints = 0,
            model_size = 0,
            solver_name = "Gurobi",
            solver_iterations = 0,
            solver_nodes = 0
            )
        filename = "Random_plan_$(i)"
        write_solution(random_plan_not_solved,foldername,filename,HPC_folder_save)
        # Punishing shifts
        mo = create_random_stowageplan_model(cs_random, not_stowed, cargoc, vessel, slots)
        set_time_limit_sec(mo, 60 * 15) # 5 minutes
        set_silent(mo) 
        optimize!(mo)
        sol_no_shifts = extract_solution(det_pro, mo)
        write_solution(sol_no_shifts, foldername, "Solution_no_shifts_$(i)", HPC_folder_save)
        #cs_test2[j,i] = cs_random
        #ballast_water_no_shifts[j,i] = sum(value.(mo[:ballast_volume]))
        #cs_test2_no_shifts[j,i] = value.(mo[:cs])
        # allows shifts
        mo = random_stowageplan_allowshifts(cs_random,cargoc,vessel,slots)
        set_time_limit_sec(mo, 60 * 15) # 5 minutes
        set_silent(mo) 
        optimize!(mo)
        sol_shifts = extract_solution(det_pro, mo)
        write_solution(sol_shifts, foldername, "Solution_shifts_$(i)", HPC_folder_save)
        #cs_test2_shifts[j,i] = value.(mo[:cs])
        #ballast_water_shifts[j,i] = sum(value.(mo[:ballast_volume]))
    end
    # Sort cargocollection - car first
    sort_order = [2, 1, 3, 4] # car, truck, machine, secu
    cargoc_carfirst = sort_cargocollection(cargoc, sort_order)
    foldername = "Random_Stowage_Plan_Cars_first"
    for i in 1:no_test
        cs_random_carfirst, not_stowaged_carfirst = random_stowage_plan(cargoc_carfirst, slots)
        random_plan_not_solved = Solution(
            gap = Inf,
            status = 0,
            objective = Inf,
            time = 0,
            cargo_weight = 0.0,
            total_weight = 0.0,
            ballast_weight = 0.0, 
            area_utilization = 0.0,
            cargo = CargoPlacement[],
            cs = cs_random_carfirst,
            n_cargo_total = length(cargoc),
            n_cargo_loaded = length(cargoc)-length(not_stowaged_carfirst),
            shear_force = Float64[],
            bending_moment = Float64[],
            ballast_volume = Float64[],
            lcg = 0.0,
            tcg = 0.0, 
            vcg = 0.0,
            n_variables = 0,
            n_binary_variables = 0,
            n_constraints = 0,
            model_size = 0,
            solver_name = "Gurobi",
            solver_iterations = 0,
            solver_nodes = 0
            )
        write_solution(random_plan_not_solved,foldername,"Random_plan_$(i)",HPC_folder_save)
        # Punishing shifts
        mo = create_random_stowageplan_model(cs_random_carfirst, not_stowaged_carfirst, cargoc, vessel, slots)
        set_time_limit_sec(mo, 60 * 15) # 5 minutes
        set_silent(mo)
        optimize!(mo)
        sol_no_shifts = extract_solution(det_pro, mo)
        write_solution(sol_no_shifts, foldername, "Solution_no_shifts_$(i)", HPC_folder_save)
        #cs_test2_carfirst[j,i] = cs_random_carfirst
        #ballast_water_no_shifts_carfirst[j,i] = sum(value.(mo[:ballast_volume]))
        #cs_test2_no_shifts_carfirst[j,i] = value.(mo[:cs])
        # allows shifts
        mo = random_stowageplan_allowshifts(cs_random_carfirst,cargoc,vessel,slots)
        set_time_limit_sec(mo, 60 * 15) # 5 minutes
        set_silent(mo)
        optimize!(mo)
        sol_shifts = extract_solution(det_pro, mo)
        write_solution(sol_shifts, foldername, "Solution_shifts_$(i)", HPC_folder_save)
        #cs_test2_shifts_carfirst[j,i] = value.(mo[:cs])
        #ballast_water_shifts_carfirst[j,i] = sum(value.(mo[:ballast_volume]))
    end
    # Sort cargocollection - Secu first
    sort_order = [4, 1, 2, 3] # secu, truck, car, machine
    cargoc_secufirst = sort_cargocollection(cargoc, sort_order)
    foldername = "Random_Stowage_Plan_Secu_first"
    for i in 1:no_test
        cs_random_secufirst, not_stowaged_secufirst = random_stowage_plan(cargoc_secufirst, slots)
        random_plan_not_solved = Solution(
            gap = Inf,
            status = 0,
            objective = Inf,
            time = 0,
            cargo_weight = 0.0,
            total_weight = 0.0,
            ballast_weight = 0.0, 
            area_utilization = 0.0,
            cargo = CargoPlacement[],
            cs = cs_random_secufirst,
            n_cargo_total = length(cargoc),
            n_cargo_loaded = length(cargoc)-length(not_stowaged_secufirst),
            shear_force = Float64[],
            bending_moment = Float64[],
            ballast_volume = Float64[],
            lcg = 0.0,
            tcg = 0.0, 
            vcg = 0.0,
            n_variables = 0,
            n_binary_variables = 0,
            n_constraints = 0,
            model_size = 0,
            solver_name = "Gurobi",
            solver_iterations = 0,
            solver_nodes = 0
            )
        filename = "Random_plan_$(i)"
        write_solution(random_plan_not_solved,foldername,filename,HPC_folder_save)
        # Punishing shifts
        mo = create_random_stowageplan_model(cs_random_secufirst, not_stowaged_secufirst, cargoc, vessel, slots)
        set_time_limit_sec(mo, 60 * 15) # 5 minutes
        set_silent(mo)
        optimize!(mo)
        sol_no_shifts = extract_solution(det_pro, mo)
        write_solution(sol_no_shifts, foldername, "Solution_no_shifts_$(i)", HPC_folder_save)
        #cs_test2_secufirst[j,i] = cs_random_secufirst
        #ballast_water_no_shifts_secufirst[j,i] = sum(value.(mo[:ballast_volume]))
        #cs_test2_no_shifts_secufirst[j,i] = value.(mo[:cs])
        # allows shifts
        mo = random_stowageplan_allowshifts(cs_random_secufirst,cargoc,vessel,slots)
        set_time_limit_sec(mo, 60 * 15) # 5 minutes
        set_silent(mo)
        optimize!(mo)
        sol_shifts = extract_solution(det_pro, mo)
        write_solution(sol_shifts, foldername, "Solution_shifts_$(i)", HPC_folder_save)
        #cs_test2_shifts_secufirst[j,i] = value.(mo[:cs])
        #ballast_water_shifts_secufirst[j,i] = sum(value.(mo[:ballast_volume]))
    end
end
#####################
# What issues does an empty ship have.
test_problem_name = Finlandia_test[1]
HPC_folder_load = HPC_folders[1]
HPC_folder_save = HPC_folder_load
# Load solution
det_sol = get_solution_deterministic("Finlandia_deterministic",
        "Deterministic_Solution", HPC_folder_load)
det_pro = load_data(problemname1, test_problem_name, problemname3)
cs_empty = zeros(length(det_pro.cargo),length(det_pro.slots))
model_empty_no_slack = second_stage_model(cs_empty, det_pro)
set_time_limit_sec(model_empty_no_slack, 60 * 60) # 60 minutes
set_silent(model_empty_no_slack)
optimize!(model_empty_no_slack)
temp_sol = Solution(
            gap = Inf,
            status = 0,
            objective = Inf,
            time = solve_time(model_empty_no_slack),
            cargo_weight = 0.0,
            total_weight = 0.0,
            ballast_weight = 0.0, 
            area_utilization = 0.0,
            cargo = CargoPlacement[],
            cs = zeros(length(det_pro.cargo), length(det_pro.slots)),
            n_cargo_total = length(det_pro.cargo),
            n_cargo_loaded = 0,
            shear_force = Float64[],
            bending_moment = Float64[],
            ballast_volume = Float64[],
            lcg = 0.0,
            tcg = 0.0, 
            vcg = 0.0,
            n_variables = num_variables(model_empty_no_slack),
            n_binary_variables = count(is_binary, all_variables(model_empty_no_slack)),
            n_constraints = num_constraints(model_empty_no_slack; count_variable_in_set_constraints=true),
            model_size = Base.summarysize(model_empty_no_slack),
            solver_name = string(JuMP.solver_name(model_empty_no_slack)),
            solver_iterations = 0,
            solver_nodes = 0
        )
temp_sol.cs
sol = get_solution_second_stage_deterministic(det_pro, model_empty_no_slack, temp_sol)
# Save Solution for empty ship
write_solution(sol,"Finlandia_deterministic","Empty_ship",HPC_folder_save)
#############
=#