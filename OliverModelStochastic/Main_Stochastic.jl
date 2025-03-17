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
problems_sto_gen = []
problems_sto_monte = []
models_sto_gen = []
models_sto_monte = []
# solve models
solutions_sto_gen = []
solutions_sto_monte = []
for i in 1:length(scenarios)
    for j in 1:length(n_cargo_unknownweight)
        # generic method
        println("gen, s = $i, n = $j")
        pro = create_stochastic_problem(problem_det, scenarios[i], n_cargo_unknownweight[j], []) 
        push!(problems_sto_gen, pro)
        mo = create_model_stochastic(pro)
        push!(models_sto_gen, mo)
        set_silent(mo) # removes terminal output
        set_time_limit_sec(mo, time_limit)
        optimize!(mo)
        push!(solutions_sto_gen, extract_stochastic_solution(pro, mo))
        # monte carlo method
        println("Monte, s = $i, n = $j")
        pro = create_stochastic_problem(problem_det, scenarios[i], n_cargo_unknownweight[j], [], Monto_Carlo_sampling) 
        push!(problems_sto_monte, pro)
        mo = create_model_stochastic(pro)
        push!(models_sto_monte, mo)
        set_silent(mo) # removes terminal output
        set_time_limit_sec(mo, time_limit)
        optimize!(mo)
        push!(solutions_sto_monte, extract_stochastic_solution(pro, mo))
    end
end

# Get solution when knowing all cargo weights

# plot solutions
# gen
diff_ballast, diff_cargo_weight, diff_placed, diff_index, same_index = compare_solutions_print(solutions_sto_gen[1], solution_det)
plot_ballast_cargo_weight_diff(diff_ballast, diff_cargo_weight,scenarios[1],n_cargo_unknownweight[1])
# monte
diff_ballast, diff_cargo_weight, diff_placed, diff_index, same_index = compare_solutions_print(solutions_sto_monte[1], solution_det)
plot_ballast_cargo_weight_diff(diff_ballast, diff_cargo_weight,scenarios[1],n_cargo_unknownweight[1])


#lm = LoadMaster("http://localhost:5000")
#v = getcurrentvessel(lm)
#resetvessel(lm)

#upload_solution(lm, solution)
