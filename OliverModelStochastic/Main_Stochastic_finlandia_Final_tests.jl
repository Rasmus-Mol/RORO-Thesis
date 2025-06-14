# Script to run on HPC
# Has noise, uses scenario reduction, and also finds EVP solution
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
#include("src/representation/ScenarioReduction.jl")
using .StowagePlannerStochastic
using JuMP

parse_index = parse(Int, ARGS[1]) # Jobindex input 
#parse_index = 1
# Choose instance:
#test_problem_name = Finlandia_test[1]
test_problem_name = Finlandia_test[parse_index]

problem_det = load_data("finlandia", test_problem_name, "hazardous")

# Should be decided more precisely at some point
initial_scenarios = 100
# Folder name for results - date and hour
HPC_folder = "Finlandia_"*test_problem_name*"_Final_Tests_init_scenarios_$(initial_scenarios)_"*Dates.format(now(), "dd_mm_HH")
 # Describe tests if necessary
extra_info = "Ship: Finlandia, Test problem: "*test_problem_name*" - Yes Scenario reduction, Yes noise, Yes EVP"

# First job index - create problem 
if parse_index == 1
    # Creates Deterministic problem and model
    model_det = create_model(problem_det)
    set_silent(model_det) # removes terminal output
    set_time_limit_sec(model_det, 60 * 60) # 1 hour
    optimize!(model_det)
    solution_det = extract_solution(problem_det, model_det)
    # Save Solution for deterministic problem
    write_solution(solution_det,"Finlandia_deterministic","Deterministic_Solution",HPC_folder)
end

# Scenarios we test - length has to fit to
# TODO: Should change when I know how many is needed
#scenarios = [10,20,30,40,50]
#scenarios = [10,20]
scenarios = [10]
sc = scenarios[1]
#sc = scenarios[parse_index] # number scenarios for current job
n_cargo_unknownweight = [length(problem_det.cargo)] # all cargo weights are unknown
time_limit = 60 * 60 # 1 hour
repetitions = 10 # number of repetitions of same inputs

# Check if folder and file has been created, otherwise create
file_check = "Results/"*HPC_folder*"HPC_data.json"
if !isfile(file_check)
    # Save scenario and number of unknown weights
    write_HPC_data(repetitions, scenarios, n_cargo_unknownweight, time_limit, HPC_folder, extra_info)
end

# checking everything works
#=
temp = add_white_noise_to_test_instance(problem_det)
temp_pro = create_stochastic_problem(temp, initial_scenarios, n_cargo_unknownweight[1], [],Bootstrap_bookedweight_quantile) 
# TODO: how many scenarios reduce to, i.e. sc
temp_pro_end = scenario_reduced(temp_pro, sc, scenario_reduction_clustering,60)
=#

# Run tests
for i in 1:repetitions
    # uniform random sampling method
    # not necessary to do scenario reduction on this method
    pro = create_stochastic_problem(problem_det, sc, n_cargo_unknownweight[1], []) 
    # Save problem
    foldername = "Stochastic_random_sampling_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)"
    write_problem_stochastic(pro,foldername,"Stochastic_Problem",HPC_folder)
    mo = create_model_stochastic(pro)
    set_silent(mo) # removes terminal output
    set_time_limit_sec(mo, time_limit)
    optimize!(mo)
    sol = extract_stochastic_solution(pro,mo)
    # Save solution
    write_solution_stochastic(sol,foldername,"Stochastic_Solution",HPC_folder)
    cs_sol = sol.cs
    # Solve second stage when we know unknown weights
    second_stage_m = second_stage_model(cs_sol, problem_det)
    set_silent(second_stage_m) # removes terminal output
    set_time_limit_sec(second_stage_m, time_limit) 
    optimize!(second_stage_m)
    fitted_sol = get_solution_second_stage_stochastic(problem_det, second_stage_m, sol)
    # Save fitted solution
    write_solution(fitted_sol,foldername,"Fitted_Solution",HPC_folder)
    # Check if second-stage problem was feasible
    #println(fitted_sol.status)
    #if fitted_sol.status == "INFEASIBLE"
    if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        set_silent(second_stage_m_slacked) # removes terminal output
        set_time_limit_sec(second_stage_m_slacked, time_limit) 
        optimize!(second_stage_m_slacked)
        fitted_sol_slacked = get_solution_second_stage_stochastic(problem_det, second_stage_m_slacked, sol)
        write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
    end
    # EVP method for uniform random sampling method 
    foldername = "EVP_random_sampling_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)"
    pro_sto = pro
    pro = expected_value_problem(pro)
    # Save problem
    write_problem(pro,foldername,"EVP_Problem",HPC_folder)
    mo = create_model(pro)
    set_silent(mo) # removes terminal output
    set_time_limit_sec(mo, time_limit)
    optimize!(mo)
    sol = extract_solution(pro, mo)
    # Save solution
    write_solution(sol,foldername,"EVP_Solution",HPC_folder)
    cs_sol = sol.cs
    # Finding solutions to each scenario for EVP
    for l in 1:sc
        # create problems
        pro_temp = StowageProblem(
            vessel = problem_det.vessel,
            slots = problem_det.slots,
            cargo = pro_sto.cargo.items[l],   
            # Problem metadata
            name = "EVP Random sampling $(l)",
            timestamp = now()
            )
        EVP_model = second_stage_model(cs_sol, pro_temp)
        set_silent(EVP_model) # removes terminal output
        set_time_limit_sec(EVP_model, time_limit)
        optimize!(EVP_model)
        fitted_sol_EVP = get_solution_second_stage_deterministic(pro_temp, EVP_model, sol)
        write_solution(fitted_sol_EVP, foldername, "Fitted_Solution_EVP_$(l)", HPC_folder)
    end
    second_stage_m = second_stage_model(cs_sol, problem_det)
    set_silent(second_stage_m) # removes terminal output
    set_time_limit_sec(second_stage_m, time_limit) 
    optimize!(second_stage_m)
    fitted_sol = get_solution_second_stage_deterministic(problem_det, second_stage_m, sol)
    # Save fitted solution
    write_solution(fitted_sol,foldername,"Fitted_Solution",HPC_folder)
    # Check if second-stage problem was feasible
    #fitted_sol.status == "INFEASIBLE"
    if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        set_silent(second_stage_m_slacked) # removes terminal output
        set_time_limit_sec(second_stage_m_slacked, time_limit) 
        optimize!(second_stage_m_slacked)
        fitted_sol_slacked = get_solution_second_stage_deterministic(problem_det, second_stage_m_slacked, sol)
        write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
    end
    # Bootstrap method 1
    # TODO: scenario reduction
    problem_det_noise = add_white_noise_to_test_instance(problem_det)
    pro = create_stochastic_problem(problem_det_noise, initial_scenarios, n_cargo_unknownweight[1], [],Bootstrap_bookedweight_quantile) 
    # Scenario reduction
    pro = scenario_reduced(pro, sc, scenario_reduction_clustering,60)
    # Save problem
    foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)"
    write_problem_stochastic(pro,foldername,"Stochastic_Problem",HPC_folder)
    mo = create_model_stochastic(pro)
    set_silent(mo) # removes terminal output
    set_time_limit_sec(mo, time_limit)
    optimize!(mo)
    sol = extract_stochastic_solution(pro,mo)
    # Save solution
    write_solution_stochastic(sol,foldername,"Stochastic_Solution",HPC_folder)
    cs_sol = sol.cs
    # Solve second stage when we know unknown weights
    second_stage_m = second_stage_model(cs_sol, problem_det)
    set_silent(second_stage_m) # removes terminal output
    set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
    optimize!(second_stage_m)
    fitted_sol = get_solution_second_stage_stochastic(problem_det, second_stage_m, sol)
    # Save fitted solution
    write_solution(fitted_sol,foldername,"Fitted_Solution",HPC_folder)
    # Check if second-stage problem was feasible
    #if fitted_sol.status == "INFEASIBLE"
    if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        set_silent(second_stage_m_slacked) # removes terminal output
        set_time_limit_sec(second_stage_m_slacked, time_limit) 
        optimize!(second_stage_m_slacked)
        fitted_sol_slacked = get_solution_second_stage_stochastic(problem_det, second_stage_m_slacked, sol)
        write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
    end
    # EVP method for Bootstrap method 1
    foldername = "EVP_Bootstrap1_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)"
    pro_sto = pro
    pro = expected_value_problem(pro)
    # Save problem
    write_problem(pro,foldername,"EVP_Problem",HPC_folder)
    mo = create_model(pro)
    set_silent(mo) # removes terminal output
    set_time_limit_sec(mo, time_limit)
    optimize!(mo)
    sol = extract_solution(pro, mo)
    # Save solution
    write_solution(sol,foldername,"EVP_Solution",HPC_folder)
    cs_sol = sol.cs
    # Finding solutions to each scenario for EVP
    for l in 1:sc
        # create problems
        pro_temp = StowageProblem(
            vessel = problem_det.vessel,
            slots = problem_det.slots,
            cargo = pro_sto.cargo.items[l],   
            # Problem metadata
            name = "EVP Bootstrap sampling $(l)",
            timestamp = now()
            )
        EVP_model = second_stage_model(cs_sol, pro_temp)
        set_silent(EVP_model) # removes terminal output
        set_time_limit_sec(EVP_model, time_limit)
        optimize!(EVP_model)
        fitted_sol_EVP = get_solution_second_stage_deterministic(pro_temp, EVP_model, sol)
        write_solution(fitted_sol_EVP, foldername, "Fitted_Solution_EVP_$(l)", HPC_folder)
    end
    second_stage_m = second_stage_model(cs_sol, problem_det)
    set_silent(second_stage_m) # removes terminal output
    set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
    optimize!(second_stage_m)
    fitted_sol = get_solution_second_stage_deterministic(problem_det, second_stage_m, sol)
    # Save fitted solution
    write_solution(fitted_sol,foldername,"Fitted_Solution",HPC_folder)
    # Check if second-stage problem was feasible
    #fitted_sol.status == "INFEASIBLE"
    if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        set_silent(second_stage_m_slacked) # removes terminal output
        set_time_limit_sec(second_stage_m_slacked, time_limit) 
        optimize!(second_stage_m_slacked)
        fitted_sol_slacked = get_solution_second_stage_deterministic(problem_det, second_stage_m_slacked, sol)
        write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
    end
end

