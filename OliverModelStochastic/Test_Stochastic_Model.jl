# Script to test the stochastic model
push!(LOAD_PATH, pwd())
# Doens't work with out this, dont know why
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

# Creates Stochastic problem and model
# parameters
scenarios = [10,20,50]
n_cargo_unknownweight = [10,50,length(problem_det.cargo)] # all cargo weights are unknown
time_limit = 60 * 5
# create problems and models
problems_sto_gen = create_stochastic_problem(problem_det, scenarios[1], n_cargo_unknownweight[1],[])
problems_sto_monte = create_stochastic_problem(problem_det, scenarios[1], n_cargo_unknownweight[1],[],Monto_Carlo_sampling)
models_sto_gen = create_model_stochastic(problems_sto_gen)
models_sto_monte = create_model_stochastic(problems_sto_monte)
# solve models
set_silent(models_sto_gen) # removes terminal output
set_silent(models_sto_monte) # removes terminal output
set_time_limit_sec(models_sto_gen, 60 * 5) # 5 minutes to solve model
set_time_limit_sec(models_sto_monte, 60 * 5) # 5 minutes to solve model
optimize!(models_sto_gen)
optimize!(models_sto_monte)
solutions_sto_gen = extract_stochastic_solution(problems_sto_gen, models_sto_gen)
solutions_sto_monte = extract_stochastic_solution(problems_sto_monte, models_sto_monte)
cs_gen = solutions_sto_gen.cs
cs_monte = solutions_sto_monte.cs

# Get solution, when knowing all cargo weights
new_model_gen = second_stage_model(cs_gen, problem_det)
new_model_monte = second_stage_model(cs_monte, problem_det)
set_silent(new_model_gen) # removes terminal output
set_silent(new_model_monte) # removes terminal output
set_time_limit_sec(new_model_gen, 60 * 5) # 5 minutes to solve model
set_time_limit_sec(new_model_monte, 60 * 5) # 5 minutes to solve model
optimize!(new_model_gen)
optimize!(new_model_monte)
fitted_solution_gen = get_solution_second_stage(problem_det, new_model_gen, solutions_sto_gen)
fitted_solution_monte = get_solution_second_stage(problem_det, new_model_monte, solutions_sto_monte)


# plot solutions and problems 
#using Plots
#include("src/plots/weight_plots.jl")
# deterministic problem
cargo_plots = plot_cargo_OG(problem_det)
#display(cargo_plots[1])
#display(cargo_plots[2])
# gen
diff_ballast, diff_cargo_weight, diff_placed, diff_index, same_index = compare_solutions_print(solutions_sto_gen, solution_det,1)
plot_ballast_weight_diff(diff_ballast,scenarios[1],n_cargo_unknownweight[1])
gen_scenarios_plots = plot_cargo_weights(problems_sto_gen,[i for i in 1:min(scenarios[1],10)])
#display(gen_scenarios_plots[1])
# monte
diff_ballast, diff_cargo_weight, diff_placed, diff_index, same_index = compare_solutions_print(solutions_sto_monte, solution_det,1)
plot_ballast_weight_diff(diff_ballast,scenarios[1],n_cargo_unknownweight[1])
monte_scenarios_plots = plot_cargo_weights(problems_sto_monte,[i for i in 1:min(scenarios[1],10)])
#display(monte_scenarios_plots[1])


#lm = LoadMaster("http://localhost:5000")
#v = getcurrentvessel(lm)
#resetvessel(lm)

#upload_solution(lm, solution)

