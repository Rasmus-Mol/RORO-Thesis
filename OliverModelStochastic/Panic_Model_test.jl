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
HPC_folder = "Finlandia_temp_model_all_instances_2"#_test_instance_$(parse_index)_"*Dates.format(now(), "dd_mm_HH")

 # Creates Deterministic problem and model - old slope coefficients
 model_det_1 = create_model_temp(problem_det)
 set_silent(model_det_1) # removes terminal output
 set_time_limit_sec(model_det_1, 60 * 60) # 1 hour
 optimize!(model_det_1)
 solution_det = extract_solution(problem_det, model_det_1)
 # Save Solution for deterministic problem
 write_solution(solution_det,test_problem_name,"Deterministic_Solution_old_slope",HPC_folder)

 using JSON, JSON3
temp_vcv_ballast = collect(value.(model_det_1[:vcg_ballast]))
temp_vcg_total = collect(value.(model_det_1[:vcg_total]))
temp_zmin = collect(value.(model_det_1[:z_min]))
z_min = findall(x -> x == 1, temp_zmin)
temp_zmax = collect(value.(model_det_1[:z_max]))
z_max = findall(x -> x == 1, temp_zmax)
temp = [temp_vcv_ballast,temp_vcg_total, z_min, z_max]
open(joinpath("Results",HPC_folder,test_problem_name,"VCGs_old_slope.json"), "w") do file
	JSON.print(file, temp, 4)  # Pretty-print with indentation
end

  # Creates Deterministic problem and model - new wrong slope coefficients
  model_det_2 = create_model_test(problem_det)
  set_silent(model_det_2) # removes terminal output
  set_time_limit_sec(model_det_2, 60 * 60) # 1 hour
  optimize!(model_det_2)
  solution_det = extract_solution(problem_det, model_det_2)
  # Save Solution for deterministic problem
  write_solution(solution_det,test_problem_name,"Deterministic_Solution_my_slope",HPC_folder)
  temp1_vcg_ballast = collect(value.(model_det_2[:vcg_ballast]))
  temp1_vcg_total = collect(value.(model_det_2[:vcg_total]))
  temp1_zmin = collect(value.(model_det_2[:z_min]))
  z_min1 = findall(x -> x == 1, temp1_zmin)
  temp1_zmax = collect(value.(model_det_2[:z_max]))
  z_max1 = findall(x -> x == 1, temp1_zmax)
  temp1 = [temp1_vcg_ballast,temp1_vcg_total,z_min1,z_max1]
  open(joinpath("Results",HPC_folder,test_problem_name,"VCGs_my_slope.json"), "w") do file
	  JSON.print(file, temp1, 4)  # Pretty-print with indentation
  end

slope_my_wrong = calculate_vcg_slopes(problem_det.vessel)
slope_old_inital_slope = calculate_vcg_slopes_temp(problem_det.vessel)

open(joinpath("Results",HPC_folder,test_problem_name,"old_slope.json"), "w") do file
	JSON.print(file, slope_old_inital_slope, 4)  # Pretty-print with indentation
end
open(joinpath("Results",HPC_folder,test_problem_name,"my_slope.json"), "w") do file
	JSON.print(file, slope_my_wrong, 4)  # Pretty-print with indentation
end