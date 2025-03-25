# Script to run on HPC
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
using .StowagePlannerStochastic
using JuMP


# Creates Deterministic problem and model
problem_det = load_data("finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous")
model_det = create_model(problem_det)
set_silent(model_det) # removes terminal output
set_time_limit_sec(model_det, 60 * 5)
optimize!(model_det)
solution_det = extract_solution(problem_det, model_det)
# Save Solution
write_solution(solution_det,"Finlandia_deterministic","Deterministic_Solution")



# Creates Stochastic problem and model
# parameters
#scenarios = [10,20,50]
scenarios = [10]
sc = length(scenarios)
#n_cargo_unknownweight = [10,50,length(problem_det.cargo)] # all cargo weights are unknown
n_cargo_unknownweight = [10]
n = length(n_cargo_unknownweight)
time_limit = 60 * 5
repetitions = 1 # number of repetitions of same inputs
# create problems and models
problems_sto_gen = Array{Any}(undef, repetitions, sc,n)
problems_EVP = Array{Any}(undef, repetitions, sc,n)
models_sto_gen = Array{Any}(undef, repetitions, sc,n)
models_EVP = Array{Any}(undef, repetitions, sc,n)
# solve models
solutions_sto_gen = Array{Any}(undef, repetitions, sc,n)
solutions_EVP = Array{Any}(undef, repetitions, sc,n)
cs_gen = Array{Any}(undef, repetitions, sc,n)
cs_EVP = Array{Any}(undef, repetitions, sc,n)
fitted_sol_gen = Array{Any}(undef, repetitions, sc,n)
fitted_sol_EVP = Array{Any}(undef, repetitions, sc,n)

for i in 1:repetitions
    println("Iteration: ", i)
    for j in 1:sc
        for k in 1:n
            # generic method
            pro = create_stochastic_problem(problem_det, scenarios[j], n_cargo_unknownweight[k], []) 
            problems_sto_gen[i,j,k] = pro
            mo = create_model_stochastic(pro)
            models_sto_gen[i,j,k] = mo # Why save the model?
            set_silent(mo) # removes terminal output
            set_time_limit_sec(mo, time_limit)
            optimize!(mo)
            sol = extract_stochastic_solution(pro,mo)
            solutions_sto_gen[i,j,k] = sol
            cs_sol = sol.cs
            cs_gen[i,j,k] = cs_sol
            second_stage_m = second_stage_model(cs_sol, problem_det)
            set_silent(second_stage_m) # removes terminal output
            set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
            optimize!(second_stage_m)
            fitted_sol = get_solution_second_stage(problem_det, second_stage_m, sol)
            fitted_sol_gen[i,j,k] = fitted_sol

            # EVP method
            pro = expected_value_problem(pro)
            problems_EVP[i,j,k] = pro
            mo = create_model(pro)
            models_EVP[i,j,k] = mo # Why save the model?
            set_silent(mo) # removes terminal output
            set_time_limit_sec(mo, time_limit)
            optimize!(mo)
            sol = extract_solution(pro, mo)
            solutions_sto_gen[i,j,k] = sol
            cs_sol = sol.cs
            cs_EVP[i,j,k] = cs_sol
            second_stage_m = second_stage_model(cs_sol, problem_det)
            set_silent(second_stage_m) # removes terminal output
            set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
            optimize!(second_stage_m)
            fitted_sol = get_solution_second_stage(problem_det, second_stage_m, sol)
            fitted_sol_EVP[i,j,k] = fitted_sol
        end
    end
end

