# Main script to test if we slack deck weight limit what then happens
# TODO
# Script to run on HPC
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
#include("src/utils/addnoise.jl")
using .StowagePlannerStochastic
using JuMP

parse_index = parse(Int, ARGS[1]) # Jobindex input
#parse_index = 1
# Choose instance:
test_problem_name = Finlandia_test[8]

problem_det = load_data("finlandia", test_problem_name, "hazardous")

# Folder name for results - date and hour
HPC_folder = "Finlandia_"*test_problem_name*"_"*Dates.format(now(), "dd_mm_HH")
 # Describe tests if necessary
extra_info = "Ship: Finlandia, Test problem: "*test_problem_name*" - No Scenario reduction, Yes noise, No EVP"

# First job index - create problem 
# TODO:
# Should also be slacked
if parse_index == 1
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
repetitions = 10 # number of repetitions of same inputs

# Check if folder and file has been created, otherwise create
file_check = "Results/"*HPC_folder*"HPC_data.json"
if !isfile(file_check)
    # Save scenario and number of unknown weights
    write_HPC_data(repetitions, scenarios, n_cargo_unknownweight, time_limit, HPC_folder, extra_info)
end