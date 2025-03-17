push!(LOAD_PATH, pwd())
using StowagePlanner
using JuMP

problem = load_data("finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous")
model = create_model(problem)

set_silent(model) # removes terminal output
set_time_limit_sec(model, 60 * 5)
optimize!(model)
solution = extract_solution(problem, model)

plot_solution(solution)

#lm = LoadMaster("http://localhost:5000")
#v = getcurrentvessel(lm)
#resetvessel(lm)

#upload_solution(lm, solution)
