# Script to run on HPC
#=
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
include("src/utils/SaveData.jl")
#include("src/utils/addnoise.jl")
using .StowagePlannerStochastic
using JuMP
=#
include("packages_and_files.jl")

HPC_folder_load = "Finlandia_no_cars_heavy_100_28_05_19"
test_instance = Finlandia_test[8]
problem_det = load_data("finlandia", test_instance, "hazardous")
# load problem
problemname1, problemname2, problemname3 = "finlandia", test_instance, "hazardous"
repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder_load)
foldername = "Stochastic_Bootstrap1_rep$(9)_sc$(scenarios[5])_unknown$(n_unknown[1])_time$(time_limit)"
filename = "Stochastic_Problem"
pro = get_stochastic_problem(foldername,
filename,HPC_folder_load,problemname1,problemname2,problemname3)

# Try a more time to see if it gets better
time_limit = 60*60*5 # 5 hours
foldername = "Stochastic_Bootstrap1_rep$(9)_sc$(scenarios[5])_unknown$(n_unknown[1])_time$(time_limit)"


HPC_folder = "Small_test_temp"
mo = create_model_stochastic(pro)
set_silent(mo) # removes terminal output
set_time_limit_sec(mo, time_limit)
optimize!(mo)
# Save solution
sol = extract_stochastic_solution(pro,mo)
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
# Check if second-stage problem was feasible
#println(fitted_sol.status)
#if fitted_sol.status == "INFEASIBLE"
if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
    second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
    set_silent(second_stage_m_slacked) # removes terminal output
    set_time_limit_sec(second_stage_m_slacked, time_limit) 
    optimize!(second_stage_m_slacked)
    fitted_sol_slacked = get_solution_second_stage_stochastic(problem_det, second_stage_m_slacked, sol)
    write_solution(fitted_sol_slacked,foldername,"Fitted_Solution_slacked",HPC_folder)
    write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
end
