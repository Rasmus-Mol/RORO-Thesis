# Script to run on HPC
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
using .StowagePlannerStochastic
using JuMP

parse_index = parse(Int, ARGS[1]) # Jobindex input


problem_det = load_data("finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous")

# First job index - create problem 
if parse_index == 1
    # Folder name for results - date and hour
    HPC_folder = "Finlandia_"*Dates.format(now(), "dd_mm_HH")
    # Creates Deterministic problem and model
    model_det = create_model(problem_det)
    set_silent(model_det) # removes terminal output
    set_time_limit_sec(model_det, 60 * 60) # 1 hour
    optimize!(model_det)
    solution_det = extract_solution(problem_det, model_det)
    # Save Solution for deterministic problem
    write_solution(solution_det,"Finlandia_deterministic","Deterministic_Solution",HPC_folder)
end

# Creates Stochastic problem and model
# parameters

# Scenarios we test - length has to fit to
scenarios = [10,20,30,40,50]
sc = scenarios[parse_index] # number scenarios for current job
n_cargo_unknownweight = [length(problem_det.cargo)] # all cargo weights are unknown
time_limit = 60 * 60 # 1 hour
repetitions = 1 # number of repetitions of same inputs

if parse_index == 1
    # Save scenario and number of unknown weights
    write_HPC_data(repetitions, scenarios, n_cargo_unknownweight, time_limit, HPC_folder)
end
# Run tests
for i in 1:repetitions
    # uniform random sampling method
    pro = create_stochastic_problem(problem_det, sc, n_cargo_unknownweight[1], []) 
    # Save problem
    foldername = "Stochastic_random_sampling_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[k])_time$(time_limit)"
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

    # Bootstrap method 1
    pro = create_stochastic_problem(problem_det, scenarios[j], n_cargo_unknownweight[k], [],Bootstrap_bookedweight_quantile) 
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
end


