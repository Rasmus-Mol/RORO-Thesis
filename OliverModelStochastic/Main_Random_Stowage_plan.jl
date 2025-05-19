# Script to run on HPC - Tests Random stowage plan
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
#include("src/StowagePlannerStochastic.jl")
#using .StowagePlannerStochastic
#using JuMP

include("packages_and_files.jl")



test_problem_name_cars_60 = Finlandia_test[1]
test_problem_name_cars_100 = Finlandia_test[2]
problem_det_cars_60 = load_data("finlandia", test_problem_name_cars_60, "hazardous")
problem_det_cars_100 = load_data("finlandia", test_problem_name_cars_100, "hazardous")
CargoC_cars_60 = problem_det_cars_60.cargo
CargoC_cars_100 = problem_det_cars_100.cargo
ntypes_cars_60 = [length(findall(x->x.cargo_type_id == i, CargoC_cars_60)) for i in 1:4]
ntypes_cars_100 = [length(findall(x->x.cargo_type_id == i, CargoC_cars_100)) for i in 1:4]
slots_fin = problem_det_cars_60.slots
vessel_fin = problem_det_cars_60.vessel
sort_order = [1,4,3,2]
# Random plan with cargo from 
cs_cars_60, not_stowaged_cars_60 = random_stowage_plan(CargoC_cars_60, slots_fin)
cs_cars_100, not_stowaged_cars_100 = random_stowage_plan(CargoC_cars_100, slots_fin)
println("Number of cargo not placed - cars60: ", length(not_stowaged_cars_60))
println("Number of cargo not placed - cars100: ", length(not_stowaged_cars_100))
println("Weight not placed - cars60: ", length(not_stowaged_cars_60)>0 ? sum([not_stowaged_cars_60[i].weight for i in 1:length(not_stowaged_cars_60)]) : 0)
println("Weight not placed - cars100: ", length(not_stowaged_cars_100)>0 ? sum([not_stowaged_cars_100[i].weight for i in 1:length(not_stowaged_cars_100)]) : 0)
plot_solution_random_plan(cs_cars_60, CargoC_cars_60, slots_fin)
plot_solution_random_plan(cs_cars_100, CargoC_cars_100, slots_fin)
# Sort cargocollection, so we do not start with placing cars
CargoC_cars_60_sorted = sort_cargocollection(CargoC_cars_60, sort_order)
CargoC_cars_100_sorted = sort_cargocollection(CargoC_cars_100, sort_order)
cs_cars_60_sorted, not_stowaged_cars_60_sorted = random_stowage_plan(CargoC_cars_60_sorted, slots_fin)
cs_cars_100_sorted, not_stowaged_cars_100_sorted = random_stowage_plan(CargoC_cars_100_sorted, slots_fin)
println("Number of cargo not placed, sorted - cars60: ", length(not_stowaged_cars_60_sorted))
println("Number of cargo not placed, sorted - cars100: ", length(not_stowaged_cars_100_sorted))
println("Weight not placed, sorted - cars60: ", length(not_stowaged_cars_60_sorted)>0 ? sum([not_stowaged_cars_60_sorted[i].weight for i in 1:length(not_stowaged_cars_60_sorted)]) : 0)
println("Weight not placed, sorted - cars100: ", length(not_stowaged_cars_100_sorted)>0 ? sum([not_stowaged_cars_100_sorted[i].weight for i in 1:length(not_stowaged_cars_100_sorted)]) : 0)
plot_solution_random_plan(cs_cars_60_sorted, CargoC_cars_60_sorted, slots_fin)
plot_solution_random_plan(cs_cars_100_sorted, CargoC_cars_100_sorted, slots_fin)
# Create random cargo collection
# random_cargocollection(ntypes::Vector{Int}, shuffle::Bool = true)

# Solving the problem with the random stowage plan
# Making it feasible
model_cars_60 = create_random_stowageplan_model(cs_cars_60, not_stowaged_cars_60, CargoC_cars_60, vessel_fin, slots_fin)
set_time_limit_sec(model_cars_60, 60 * 5) # 5 minutes
set_silent(model_cars_60) 
optimize!(model_cars_60)
y_60 = value.(model_cars_60[:y])
println("Number of changes:", sum(y_60))
model_cars_100 = create_random_stowageplan_model(cs_cars_100, not_stowaged_cars_100, CargoC_cars_100, vessel_fin, slots_fin)
set_time_limit_sec(model_cars_100, 60 * 5) # 5 minutes
set_silent(model_cars_100) 
optimize!(model_cars_100)
y_100 = value.(model_cars_100[:y])
println("Number of changes:", sum(y_100))
model_cars_60_sorted = create_random_stowageplan_model(cs_cars_60_sorted, not_stowaged_cars_60_sorted, CargoC_cars_60_sorted, vessel_fin, slots_fin)
set_time_limit_sec(model_cars_60_sorted, 60 * 5) # 5 minutes
set_silent(model_cars_60_sorted) 
optimize!(model_cars_60_sorted)
y_60_sorted = value.(model_cars_60_sorted[:y])
println("Number of changes:", sum(y_60_sorted))
model_cars_100_sorted = create_random_stowageplan_model(cs_cars_100_sorted, not_stowaged_cars_100_sorted, CargoC_cars_100_sorted, vessel_fin, slots_fin)
set_time_limit_sec(model_cars_100_sorted, 60 * 5) # 5 minutes
set_silent(model_cars_100_sorted) 
optimize!(model_cars_100_sorted)
y_100_sorted = value.(model_cars_100_sorted[:y])
println("Number of changes:", sum(y_100_sorted))


# create_random_stowageplan_model(cs_old,not_stowed,cargo,vessel,slots, new_cargo_allowed::Bool = false)