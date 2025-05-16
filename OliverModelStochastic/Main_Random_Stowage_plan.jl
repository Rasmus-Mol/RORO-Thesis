# Script to run on HPC - Tests Random stowage plan
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
using .StowagePlannerStochastic
using JuMP


test_problem_name = Finlandia_test[1]
problem_det = load_data("finlandia", test_problem_name, "hazardous")
CargoC = problem_det.cargo
slots = problem_det.slots
sort_order = [1,4,3,2]
cs, not_stowaged = random_stowage_plan(CargoC, slots)


random_cargocollection(ntypes::Vector{Int}, shuffle::Bool = true)
create_random_stowageplan_model(cs_old,not_stowed,cargo,vessel,slots, new_cargo_allowed::Bool = false)