# Script to run on HPC
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
using .StowagePlannerStochastic
using JuMP

# Folder name for results - date and time
HPC_folder = "Finlandia_"*Dates.format(now(), "dd_mm_HH_MM_SS")

# Creates Deterministic problem and model
problem_det = load_data("finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous")
model_det = create_model(problem_det)
set_silent(model_det) # removes terminal output
set_time_limit_sec(model_det, 60 * 5)
optimize!(model_det)
solution_det = extract_solution(problem_det, model_det)
# Save Solution for deterministic problem
write_solution(solution_det,"Finlandia_deterministic","Deterministic_Solution",HPC_folder)

# Creates Stochastic problem and model
# parameters
scenarios = [10,20,30,40,50]
#scenarios = [10]
sc = length(scenarios)
n_cargo_unknownweight = [30,80,length(problem_det.cargo)] # all cargo weights are unknown
#n_cargo_unknownweight = [10]
n = length(n_cargo_unknownweight)
time_limit = 60 * 15 # 5 minutes
repetitions = 1 # number of repetitions of same inputs

# Problems and models - Can be deleted now, no longer necessary
#problems_sto_gen = Array{Any}(undef, repetitions, sc,n)
#problems_EVP = Array{Any}(undef, repetitions, sc,n)
#models_sto_gen = Array{Any}(undef, repetitions, sc,n)
#models_EVP = Array{Any}(undef, repetitions, sc,n)
# Solution
#solutions_sto_gen = Array{Any}(undef, repetitions, sc,n)
#solutions_EVP = Array{Any}(undef, repetitions, sc,n)
#fitted_sol_gen = Array{Any}(undef, repetitions, sc,n)
#fitted_sol_EVP = Array{Any}(undef, repetitions, sc,n)

# Save scenario and number of unknown weights
write_HPC_data(repetitions, scenarios, n_cargo_unknownweight, time_limit, HPC_folder)
# Run tests
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # uniform random sampling method
            pro = create_stochastic_problem(problem_det, scenarios[j], n_cargo_unknownweight[k], []) 
            # Save problem
            foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_cargo_unknownweight[k])_time$(time_limit)"
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
            # EVP method for uniform random sampling method 
            foldername = "EVP_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_cargo_unknownweight[k])_time$(time_limit)"
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

            # Bootstrap method 1
            pro = create_stochastic_problem(problem_det, scenarios[j], n_cargo_unknownweight[k], [],Bootstrap_bookedweight_quantile) 
            # Save problem
            foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_cargo_unknownweight[k])_time$(time_limit)"
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
            foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_cargo_unknownweight[k])_time$(time_limit)"
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
    end
end


