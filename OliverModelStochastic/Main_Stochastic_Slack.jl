# Script to run on HPC
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
using .StowagePlannerStochastic
using JuMP

parse_index = parse(Int, ARGS[1]) # Jobindex input
#parse_index = 1
# Choose instance:
test_problem_name = Finlandia_test[8]

problem_det = load_data("finlandia", test_problem_name, "hazardous")
problemname1, problemname2, problemname3 = "finlandia", test_problem_name, "hazardous"
# Folder name for results - date and hour
HPC_folder_save = "Finlandia_"*test_problem_name*"_Slacktest_"*Dates.format(now(), "dd_mm_HH")
HPC_folder_load = "Finlandia_"*test_problem_name*"_15_05_09"
 # Describe tests if necessary
extra_info = "Ship: Finlandia, Test problem: "*test_problem_name*" - Slack test"

# First job index - create problem 
if parse_index == 1
    # Creates Deterministic problem and model
    model_det = create_model(problem_det)
    set_silent(model_det) # removes terminal output
    set_time_limit_sec(model_det, 60 * 60) # 1 hour
    optimize!(model_det)
    solution_det = extract_solution(problem_det, model_det)
    # Save Solution for deterministic problem
    write_solution(solution_det,"Finlandia_deterministic","Deterministic_Solution",HPC_folder_save)
end

# Creates Stochastic problem and model
# parameters

# Scenarios we test - length has to fit to
scenarios = [10,20,30,40,50]
sc = scenarios[parse_index] # number scenarios for current job
n_cargo_unknownweight = [length(problem_det.cargo)] # all cargo weights are unknown
fractions = [1,0.9, 0.8]
fraction = fractions[3] # change to test different fractions
time_limit = 60 * 60 # 1 hour
#repetitions = 5 # number of repetitions of same inputs

# load data from test
repetitions, a1, a2, a3, note = get_HPC_data(HPC_folder_load)
# Run tests
for i in 1:repetitions
    # uniform random sampling method
    # Load problem 
    #pro = create_stochastic_problem(problem_det, sc, n_cargo_unknownweight[1], [])
    foldername = "Stochastic_random_sampling_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)"
    pro = get_stochastic_problem(foldername,"Stochastic_Problem",HPC_folder_load,problemname1,problemname2,problemname3)
    # Model
    mo = create_model_stochastic_cargo_fraction(pro,fraction)
    set_silent(mo) # removes terminal output
    set_time_limit_sec(mo, time_limit)
    optimize!(mo)
    sol = extract_stochastic_solution(pro,mo)
    # Save solution
    foldername = "Stochastic_random_sampling_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)_slacktest_$(fraction)"
    write_solution_stochastic(sol,foldername,"Stochastic_Solution",HPC_folder_save)
    write_slack(HPC_folder_save, foldername, "Slack_test", mo)
    #=
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
    if fitted_sol.status == "INFEASIBLE"
        second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        set_silent(second_stage_m_slacked) # removes terminal output
        set_time_limit_sec(second_stage_m_slacked, time_limit) 
        optimize!(second_stage_m_slacked)
        fitted_sol_slacked = get_solution_second_stage_stochastic(problem_det, second_stage_m_slacked, sol)
        write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
    end
    =#
    #=
    # EVP method for uniform random sampling method 
    foldername = "EVP_random_sampling_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)"
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
    second_stage_m = second_stage_model(cs_sol, problem_det)
    set_silent(second_stage_m) # removes terminal output
    set_time_limit_sec(second_stage_m, time_limit) 
    optimize!(second_stage_m)
    fitted_sol = get_solution_second_stage_deterministic(problem_det, second_stage_m, sol)
    # Save fitted solution
    write_solution(fitted_sol,foldername,"Fitted_Solution",HPC_folder)
    # Check if second-stage problem was feasible
    if fitted_sol.status == "INFEASIBLE"
        second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        set_silent(second_stage_m_slacked) # removes terminal output
        set_time_limit_sec(second_stage_m_slacked, time_limit) 
        optimize!(second_stage_m_slacked)
        fitted_sol_slacked = get_solution_second_stage_deterministic(problem_det, second_stage_m_slacked, sol)
        write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
    end
    =#
    # Bootstrap method 1
    #pro = create_stochastic_problem(problem_det, sc, n_cargo_unknownweight[1], [],Bootstrap_bookedweight_quantile) 
    # load problem
    foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)"
    pro = get_stochastic_problem(foldername,"Stochastic_Problem",HPC_folder_load,problemname1,problemname2,problemname3)
    # Model
    mo = create_model_stochastic_cargo_fraction(pro,fraction)
    set_silent(mo) # removes terminal output
    set_time_limit_sec(mo, time_limit)
    optimize!(mo)
    sol = extract_stochastic_solution(pro,mo)
    # Save solution
    foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)_slacktest_$(fraction)"
    write_solution_stochastic(sol,foldername,"Stochastic_Solution",HPC_folder_save)
    write_slack(HPC_folder_save, foldername, "Slack_test", mo)
    #=
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
    if fitted_sol.status == "INFEASIBLE"
        second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        set_silent(second_stage_m_slacked) # removes terminal output
        set_time_limit_sec(second_stage_m_slacked, time_limit) 
        optimize!(second_stage_m_slacked)
        fitted_sol_slacked = get_solution_second_stage_stochastic(problem_det, second_stage_m_slacked, sol)
        write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
    end
    =#
    #=
    # EVP method for Bootstrap method 1
    foldername = "EVP_Bootstrap1_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)"
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
    second_stage_m = second_stage_model(cs_sol, problem_det)
    set_silent(second_stage_m) # removes terminal output
    set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
    optimize!(second_stage_m)
    fitted_sol = get_solution_second_stage_deterministic(problem_det, second_stage_m, sol)
    # Save fitted solution
    write_solution(fitted_sol,foldername,"Fitted_Solution",HPC_folder)
    # Check if second-stage problem was feasible
    if fitted_sol.status == "INFEASIBLE"
        second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        set_silent(second_stage_m_slacked) # removes terminal output
        set_time_limit_sec(second_stage_m_slacked, time_limit) 
        optimize!(second_stage_m_slacked)
        fitted_sol_slacked = get_solution_second_stage_deterministic(problem_det, second_stage_m_slacked, sol)
        write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
    end
    =#
end


