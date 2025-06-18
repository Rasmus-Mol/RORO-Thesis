push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
#include("src/utils/addnoise.jl")
using .StowagePlannerStochastic
using JuMP


parse_index = parse(Int, ARGS[1]) # Jobindex input
#parse_index = 7
# Choose instance:
test_problem_name = Finlandia_test[parse_index]

problem_det = load_data("finlandia", test_problem_name, "hazardous")
HPC_folder = "Finlandia_temp_model_all_instances"#_test_instance_$(parse_index)_"*Dates.format(now(), "dd_mm_HH")

 # Creates Deterministic problem and model - old slope coefficients
 model_det_1 = create_model_temp(problem_det)
 set_silent(model_det_1) # removes terminal output
 set_time_limit_sec(model_det_1, 60 * 60 * 3) # 2 hour
 optimize!(model_det_1)
 solution_det = extract_solution(problem_det, model_det_1)
 # Save Solution for deterministic problem
 write_solution(solution_det,"Finlandia_deterministic_old_slope"*test_problem_name,"Deterministic_Solution",HPC_folder)

 using JSON, JSON3
temp = collect(value.(model_det_1[:vcg_ballast]))
open(joinpath("Results",HPC_folder,"Finlandia_deterministic_old_slope"*test_problem_name,"VCG_ballast.json"), "w") do file
	JSON.print(file, temp, 4)  # Pretty-print with indentation
end

  # Creates Deterministic problem and model - new wrong slope coefficients
  model_det_2 = create_model_test(problem_det)
  set_silent(model_det_2) # removes terminal output
  set_time_limit_sec(model_det_2, 60 * 60 * 3) # 3 hour
  optimize!(model_det_2)
  solution_det = extract_solution(problem_det, model_det_2)
  # Save Solution for deterministic problem
  write_solution(solution_det,"Finlandia_deterministic_my_slope"*test_problem_name,"Deterministic_Solution",HPC_folder)
  temp1 = collect(value.(model_det_2[:vcg_ballast]))
  open(joinpath("Results",HPC_folder,"Finlandia_deterministic_my_slope"*test_problem_name,"VCG_ballast.json"), "w") do file
	  JSON.print(file, temp1, 4)  # Pretty-print with indentation
  end

slope_my_wrong = calculate_vcg_slopes(problem_det.vessel)
slope_old_inital_slope = calculate_vcg_slopes_temp(problem_det.vessel)

open(joinpath("Results",HPC_folder,"Finlandia_deterministic_old_slope"*test_problem_name,"old_slope.json"), "w") do file
	JSON.print(file, slope_old_inital_slope, 4)  # Pretty-print with indentation
end
open(joinpath("Results",HPC_folder,"Finlandia_deterministic_my_slope"*test_problem_name,"my_slope.json"), "w") do file
	JSON.print(file, slope_my_wrong, 4)  # Pretty-print with indentation
end