# Main script to find EVP for stochastic solutions
# Script to run on HPC
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
#include("src/utils/addnoise.jl")
using .StowagePlannerStochastic
using JuMP
HPC_folders = [
    "Finlandia_mixed_light_60_28_05_13",
    "Finlandia_mixed_light_100_28_05_13",
    "Finlandia_mixed_heavy_60_28_05_14",
    "Finlandia_mixed_heavy_100_28_05_15",
    "Finlandia_no_cars_light_60_28_05_16",
    "Finlandia_no_cars_light_100_28_05_16",
    "Finlandia_no_cars_heavy_60_28_05_19",
    "Finlandia_no_cars_heavy_100_28_05_19",
]

# Change if not Finlandia problem
problemname1, problemname3 = "finlandia", "hazardous"

parse_index = parse(Int, ARGS[1])
test_instance = Finlandia_test[parse_index]
HPC_folder = HPC_folders[parse_index]
repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
println("Extra info about test: ", note)
sc = length(scenarios)
n = length(n_unknown)
Deterministic_Solution = get_solution_deterministic("Finlandia_deterministic",
    "Deterministic_Solution", HPC_folder)
problem_det = load_data(problemname1, test_instance, problemname3)

# TODO: FIX this so i save these values
sol_evp = Array{Any}(nothing, repetitions,sc,n,scenarios[end])
obj_val_EVP = Array{Any}(nothing, repetitions,sc)
cargo_load_EVP = Array{Any}(nothing, repetitions,sc)
ballast_used_EVP = Array{Any}(nothing, repetitions,sc)
inf_index_EVP = zeros(repetitions, sc)
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            # Stochastic problem - gen
            foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            # problem
            Sto_pro_gen = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, test_instance, problemname3)
            #   EVP   
            foldername = "EVP_random_sampling_rep$(i)_sc$(sc)_unknown$(n_unknown[k])_time$(time_limit)"
            # Save stochastic problems scenarios
            pro_sto = Sto_pro_gen
            pro = expected_value_problem(pro_sto)
            # Save problem
            write_problem(pro, foldername, "EVP_Problem", HPC_folder)
            mo = create_model(pro)
            set_silent(mo) # removes terminal output
            set_time_limit_sec(mo, time_limit)
            optimize!(mo)
            sol = extract_solution(pro, mo)
            # Save solution
            write_solution(sol, foldername, "EVP_Solution", HPC_folder)
            cs_sol = sol.cs
            # Finding solutions to each scenario for EVP
            for h in 1:scenarios[j]
                # create problems
                pro_temp = StowageProblem(
                    vessel=problem_det.vessel,
                    slots=problem_det.slots,
                    cargo=pro_sto.cargo.items[h],
                    # Problem metadata
                    name="EVP Random $(h)",
                    timestamp=now()
                )
                # Just using slack model initially
                EVP_model = second_stage_model_slack(cs_sol, pro_temp)
                set_silent(EVP_model) # removes terminal output
                set_time_limit_sec(EVP_model, time_limit)
                optimize!(EVP_model)
                fitted_sol_EVP = get_solution_second_stage_deterministic(pro_temp, EVP_model, sol)
                sol_evp[i,j,k,h] = fitted_sol_EVP.objective
                #write_solution(fitted_sol_EVP, foldername, "Fitted_Solution_EVP_$(h)", HPC_folder)
                #if primal_status(EVP_model) == MOI.FEASIBLE_POINT || termination_status(EVP_model) == MOI.OPTIMAL
                #    write_slack(HPC_folder, foldername, "Fitted_Solution_$(h)", EVP_model)
                #end
            end
            #=
            # Using EVP solution with deterministic weights
            second_stage_m = second_stage_model(cs_sol, problem_det)
            set_silent(second_stage_m) # removes terminal output
            set_time_limit_sec(second_stage_m, time_limit)
            optimize!(second_stage_m)
            # Save fitted solution
            fitted_sol = get_solution_second_stage_deterministic(problem_det, second_stage_m, sol)
            write_solution(fitted_sol, foldername, "Fitted_Solution", HPC_folder)
            # Check if second-stage problem was feasible
            #if fitted_sol.status == "INFEASIBLE"
            if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
                second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
                set_silent(second_stage_m_slacked) # removes terminal output
                set_time_limit_sec(second_stage_m_slacked, time_limit)
                optimize!(second_stage_m_slacked)
                fitted_sol_slacked = get_solution_second_stage_deterministic(problem_det, second_stage_m_slacked, sol)
                write_solution(fitted_sol_slacked, foldername, "Fitted_Solution_slacked", HPC_folder)
                write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
            end
            =#
            #=
            # Stochastic problem - Boot
            foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            filename = "Stochastic_Problem"
            # problem
            Sto_pro_boot = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, test_instance, problemname3)
            # EVP   
            foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
            # Save stochastic problems scenarios
            pro_sto = Sto_pro_boot
            pro = expected_value_problem(pro_sto)
            # Save problem
            write_problem(pro, foldername, "EVP_Problem", HPC_folder)
            mo = create_model(pro)
            set_silent(mo) # removes terminal output
            set_time_limit_sec(mo, time_limit)
            optimize!(mo)
            sol = extract_solution(pro, mo)
            # Save solution
            write_solution(sol, foldername, "EVP_Solution", HPC_folder)
            cs_sol = sol.cs
            # Finding solutions to each scenario for EVP
            for h in 1:scenarios[j]
                # create problems
                pro_temp = StowageProblem(
                    vessel=problem_det.vessel,
                    slots=problem_det.slots,
                    cargo=pro_sto.cargo.items[h],
                    # Problem metadata
                    name="EVP Boot $(h)",
                    timestamp=now()
                )
                # Just using slack model initially
                EVP_model = second_stage_model_slack(cs_sol, pro_temp)
                set_silent(EVP_model) # removes terminal output
                set_time_limit_sec(EVP_model, time_limit)
                optimize!(EVP_model)
                fitted_sol_EVP = get_solution_second_stage_deterministic(pro_temp, EVP_model, sol)
                write_solution(fitted_sol_EVP, foldername, "Fitted_Solution_EVP_$(h)", HPC_folder)
                if primal_status(EVP_model) == MOI.FEASIBLE_POINT || termination_status(EVP_model) == MOI.OPTIMAL
                    write_slack(HPC_folder, foldername, "Fitted_Solution_$(h)", EVP_model)
                end
            end
            # Using EVP solution with deterministic weights
            second_stage_m = second_stage_model(cs_sol, problem_det)
            set_silent(second_stage_m) # removes terminal output
            set_time_limit_sec(second_stage_m, time_limit)
            optimize!(second_stage_m)
            # Save fitted solution
            fitted_sol = get_solution_second_stage_deterministic(problem_det, second_stage_m, sol)
            write_solution(fitted_sol, foldername, "Fitted_Solution", HPC_folder)
            # Check if second-stage problem was feasible
            #if fitted_sol.status == "INFEASIBLE"
            if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
                second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
                set_silent(second_stage_m_slacked) # removes terminal output
                set_time_limit_sec(second_stage_m_slacked, time_limit)
                optimize!(second_stage_m_slacked)
                fitted_sol_slacked = get_solution_second_stage_deterministic(problem_det, second_stage_m_slacked, sol)
                write_solution(fitted_sol_slacked, foldername, "Fitted_Solution_slacked", HPC_folder)
                write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
            end
            =#
        end
    end
end

# write matrix with objective values

