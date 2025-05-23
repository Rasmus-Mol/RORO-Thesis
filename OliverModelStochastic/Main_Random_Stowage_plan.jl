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
#parse_index = 1
# Choose instance:
test_problem_name = Finlandia_test[5]

#problem_det = load_data("finlandia", test_problem_name, "hazardous")

# Folder name for results - date and hour
HPC_folder_save = "Finlandia_"*test_problem_name*"_RandomPlan_"*Dates.format(now(), "dd_mm_HH")
 # Describe tests if necessary
extra_info = "Ship: Finlandia, Test problem: "*test_problem_name*" - Random stowage plan tests"

# Scenarios we test - length has to fit to
scenarios = [10,20,30,40,50]
sc = scenarios[parse_index] # number scenarios for current job
n_cargo_unknownweight = [length(problem_det.cargo)] # all cargo weights are unknown
time_limit = 60 * 60 # 1 hour
repetitions = 1 # number of repetitions of same inputs

# Check if folder and file has been created, otherwise create
file_check = "Results/"*HPC_folder_save*"HPC_data.json"
if !isfile(file_check)
    # Save scenario and number of unknown weights
    write_HPC_data(repetitions, scenarios, n_cargo_unknownweight, time_limit, HPC_folder, extra_info)
end

# test 1 - load solution and change weight and test feasibility.
# Load deterministic 
HPC_folder_load = "Finlandia_"*test_problem_name*"_14_05_20"
det_sol = get_solution_deterministic("Finlandia_deterministic",
"Deterministic_Solution",HPC_folder_load)
problemname1, problemname3 = "finlandia", "hazardous"
det_pro = load_data(problemname1,test_problem_name,problemname3)
det_sol_cs = Deterministic_Solution.cs
sto_pro = create_stochastic_problem(det_pro, 10, length(det_pro.cargo), []) 
f = Vector{Bool}(undef, 10)
strs = Vector{String}(undef, 10)
mo = Vector{Any}(undef, 10)
for i in 1:10
    new_c = sto_pro.cargo.items[1]
    f[i], str[i], mo[i] = feasibility_check(det_sol,det_pro,new_c)
end

# test 2 - create random stowage plan and check if it is feasible.


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