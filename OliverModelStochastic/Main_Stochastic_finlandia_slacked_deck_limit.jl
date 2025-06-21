# Script to run on HPC - Slack deck limit on deck 1
# Script to run on HPC
# Has noise, uses scenario reduction, and also finds EVP solution
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
include("src/utils/test_instances.jl")
using .StowagePlannerStochastic
using JuMP

parse_index = parse(Int, ARGS[1]) # slack-fraction input
instance_index = 8 # Which instance to test
# Choose instance:
test_problem_name = Finlandia_test[instance_index]

problem_det = load_data("finlandia", test_problem_name, "hazardous")

# Folder name for results - date and hour
HPC_folder = "Finlandia_" * test_problem_name * "_SlackDeck1_" * Dates.format(now(), "dd_mm_HH")
#slack_fraction = [[1.1, 1, 1], [1.3, 1, 1], [1.5, 1, 1]]
slack_fraction = [[1.1, 1, 1], [1.2,1,1], [1.3, 1, 1], [1.4,1,1], [1.5, 1, 1]]
# Describe tests if necessary
extra_info = "Ship: Finlandia, Test problem: " * test_problem_name * " - Yes Scenario reduction, Yes noise, Yes EVP, Yes Slack deck limit on deck: $(slack_fraction[parse_index])"

# Scenarios we test - length has to fit to
# TODO: Should change when I know how many is needed
#scenarios = [10, 20, 30, 40, 50]
scenarios = [10]#, 20]
sc = scenarios[1]#scenarios[parse_index] # number scenarios for current job
n_cargo_unknownweight = [length(problem_det.cargo)] # all cargo weights are unknown
time_limit = 60 * 60 # 1 hour
repetitions = 5 #10 # number of repetitions of same inputs

# Check if folder and file has been created, otherwise create
file_check = "Results/" * HPC_folder * "HPC_data.json"
if !isfile(file_check)
    # Save scenario and number of unknown weights
    write_HPC_data(repetitions, scenarios, n_cargo_unknownweight, time_limit, HPC_folder, extra_info)
end


# First job index - create problem 
#if parse_index == 1
    #for i in 1:length(slack_fraction)
        # Create model with slack deck limit
        model_det = create_model_slack_deck(problem_det, slack_fraction[parse_index])
        set_silent(model_det) # removes terminal output
        set_time_limit_sec(model_det, 60 * 60) # 1 hour
        optimize!(model_det)
        solution_det = extract_solution(problem_det, model_det)
        # Save Solution for deterministic problem
        write_solution(solution_det, "Finlandia_deterministic", 
        "Deterministic_Solution_slack_$(slack_fraction[parse_index])", HPC_folder)
    #end
#end

# write matrix with objective values
using JSON, JSON3
function save_EVP_objective_2D(HPC_folder::String,foldername::String,filename::String, matrix)
    if !isdir("Results/"*HPC_folder)
        mkdir("Results/"*HPC_folder)
    end
    if !isdir("Results/"*HPC_folder*"/"*foldername)
        mkdir("Results/"*HPC_folder*"/"*foldername)
    end
    # M is 2D
    M  = [[matrix[i, j] for j in 1:size(matrix, 2)] for i in 1:size(matrix, 1)]
    open(joinpath("Results",HPC_folder,foldername,filename*".json"), "w") do file
        JSON.print(file, M, 4)
    end
end 
# For later to get results
function get_obj_of_evp_2D(filename::String,foldername::String,HPC_folder::String)
    temp = open(joinpath("Results",HPC_folder,foldername,filename*".json"), "r") do file
        JSON3.read(read(file, String))#, Vector{Vector{Float64}})
    end
    rows = length(temp)
    cols = length(temp[1])
    mat = Array{Any}(undef, rows, cols)
    for i in 1:rows, j in 1:cols
        mat[i, j] = temp[i][j]
    end
    return mat
end

# Run tests
#initial_scenarios = 100
#for j in 1:length(slack_fraction)
    sol_evp_gen = Array{Any}(nothing,repetitions,sc)
    obj_val_EVP_gen = Array{Any}(nothing,repetitions,sc)
    cargo_load_EVP_gen = Array{Any}(nothing,repetitions,sc)
    ballast_used_EVP_gen = Array{Any}(nothing,repetitions,sc)
    inf_index_EVP_gen = zeros(repetitions,sc)

    sol_evp_boot = Array{Any}(nothing,repetitions,sc)
    obj_val_EVP_boot = Array{Any}(nothing,repetitions,sc)
    cargo_load_EVP_boot = Array{Any}(nothing,repetitions,sc)
    ballast_used_EVP_boot = Array{Any}(nothing,repetitions,sc)
    inf_index_EVP_boot = zeros(repetitions,sc)
    for i in 1:repetitions
        # uniform random sampling method
        # not necessary to do scenario reduction on this method
        pro = create_stochastic_problem(problem_det, sc, n_cargo_unknownweight[1], [])
        # Save problem
        foldername = "Stochastic_random_sampling_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)_slack_$(slack_fraction[parse_index])"
        # To save space
        #write_problem_stochastic(pro, foldername, "Stochastic_Problem", HPC_folder)
        #mo = create_model_stochastic(pro)
        mo = create_model_stochastic_slack_deck(pro, slack_fraction[parse_index]) 
        set_silent(mo) # removes terminal output
        set_time_limit_sec(mo, time_limit)
        optimize!(mo)
        sol = extract_stochastic_solution(pro, mo)
        # Save solution
        write_solution_stochastic(sol, foldername, "Stochastic_Solution", HPC_folder)
        cs_sol = sol.cs
        # Solve second stage when we know unknown weights
        #second_stage_m = second_stage_model(cs_sol, problem_det)
        second_stage_m = second_stage_model_deckslacked(cs_sol, problem_det, slack_fraction[parse_index])
        set_silent(second_stage_m) # removes terminal output
        set_time_limit_sec(second_stage_m, time_limit)
        optimize!(second_stage_m)
        fitted_sol = get_solution_second_stage_stochastic(problem_det, second_stage_m, sol)
        # Save fitted solution
        write_solution(fitted_sol, foldername, "Fitted_Solution", HPC_folder)
        # Check if second-stage problem was feasible
        #println(fitted_sol.status)
        #if fitted_sol.status == "INFEASIBLE"
        if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
            #second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
            second_stage_m_slacked = second_stage_model_deckslacked_slack(cs_sol, problem_det, slack_fraction[parse_index])
            set_silent(second_stage_m_slacked) # removes terminal output
            set_time_limit_sec(second_stage_m_slacked, time_limit)
            optimize!(second_stage_m_slacked)
            fitted_sol_slacked = get_solution_second_stage_stochastic(problem_det, second_stage_m_slacked, sol)
            write_solution(fitted_sol_slacked, foldername, "Fitted_Solution_slacked", HPC_folder)
            if primal_status(second_stage_m_slacked) == MOI.FEASIBLE_POINT || termination_status(second_stage_m_slacked) == MOI.OPTIMAL
                write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
            end
        end
        # EVP method for uniform random sampling method 
        foldername = "EVP_random_sampling_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)_slack_$(slack_fraction[parse_index])"
        # Save stochastic problems scenarios
        pro_sto = pro
        pro = expected_value_problem(pro)
        # Save problem
        write_problem(pro, foldername, "EVP_Problem", HPC_folder)
        #mo = create_model(pro)
        mo = create_model_slack_deck(pro, slack_fraction[parse_index])
        set_silent(mo) # removes terminal output
        set_time_limit_sec(mo, time_limit)
        optimize!(mo)
        sol = extract_solution(pro, mo)
        # Save solution
        write_solution(sol, foldername, "EVP_Solution", HPC_folder)
        cs_sol = sol.cs
        # Finding solutions to each scenario for EVP
        for l in 1:sc
            # create problems
            pro_temp = StowageProblem(
                vessel = problem_det.vessel,
                slots = problem_det.slots,
                cargo = pro_sto.cargo.items[l],   
                # Problem metadata
                name = "EVP Random $(l)",
                timestamp = now()
            )
            #EVP_model = second_stage_model_deckslacked(cs_sol, pro_temp, slack_fraction[parse_index])
            # Just create slack model and correct objevtive value later
            EVP_model = second_stage_model_deckslacked_slack(cs_sol, pro_temp, slack_fraction[parse_index])
            set_silent(EVP_model) # removes terminal output
            set_time_limit_sec(EVP_model, time_limit)
            optimize!(EVP_model)
            fitted_sol_EVP = get_solution_second_stage_deterministic(pro_temp, EVP_model, sol)
            #write_solution(fitted_sol_EVP, foldername, "Fitted_Solution_EVP_$(l)", HPC_folder)
            # Check if a feasible solution was found
            if primal_status(EVP_model) == MOI.FEASIBLE_POINT || termination_status(EVP_model) == MOI.OPTIMAL
                slack = [value.(EVP_model[:slack_deck]), value.(EVP_model[:slack_Vmax]),value.(EVP_model[:slack_Vmin]),
                    value.(EVP_model[:slack_Tmin]), value.(EVP_model[:slack_Tmax]), 
                    value.(EVP_model[:slack_Lmin]), value.(EVP_model[:slack_Lmax]),
                    #value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
                    value.(EVP_model[:slack_shearMin]),
                    value.(EVP_model[:slack_shearMax]), value.(EVP_model[:slack_bendingMax]),
                    value.(EVP_model[:slack_ballast_tanks])
                    ]
                # Check if slack is used
                if any(x->x > 0.001,vcat(slack...))
                    if sum(slack[1]) > 0.001 # deck limit
                        inf_index_EVP_gen[i,l] = 1
                    else # something else
                        inf_index_EVP_gen[i,l] = 2
                    end
                end
            else
                inf_index_EVP_gen[i,l] = -1
            end
            obj_val_EVP_gen[i,l] = fitted_sol_EVP.objective
            cargo_load_EVP_gen[i,l] = fitted_sol_EVP.n_cargo_loaded
            ballast_used_EVP_gen[i,l] = fitted_sol_EVP.ballast_weight   
        end
        # Second stage model with this solution and deterministic weights
        #second_stage_m = second_stage_model_deckslacked(cs_sol, problem_det, slack_fraction[parse_index])
        #set_silent(second_stage_m) # removes terminal output
        #set_time_limit_sec(second_stage_m, time_limit)
        #optimize!(second_stage_m)
        #fitted_sol = get_solution_second_stage_deterministic(problem_det, second_stage_m, sol)
        # Save fitted solution
        #write_solution(fitted_sol, foldername, "Fitted_Solution", HPC_folder)
        # Check if second-stage problem was feasible
        #fitted_sol.status == "INFEASIBLE"
        #if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
            #second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        #    second_stage_m_slacked = second_stage_model_deckslacked_slack(cs_sol, problem_det, slack_fraction[j])
        #    set_silent(second_stage_m_slacked) # removes terminal output
        #    set_time_limit_sec(second_stage_m_slacked, time_limit)
        #    optimize!(second_stage_m_slacked)
        #    fitted_sol_slacked = get_solution_second_stage_deterministic(problem_det, second_stage_m_slacked, sol)
        #    write_solution(fitted_sol_slacked, foldername, "Fitted_Solution_slacked", HPC_folder)
        #    if primal_status(second_stage_m_slacked) == MOI.FEASIBLE_POINT || termination_status(second_stage_m_slacked) == MOI.OPTIMAL
        #        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
        #    end 
        #end
        # Bootstrap method 1
        # TODO: scenario reduction
        problem_det_noise = add_white_noise_to_test_instance(problem_det)
        #pro = create_stochastic_problem(problem_det_noise, initial_scenarios, n_cargo_unknownweight[1], [], Bootstrap_bookedweight_quantile)
        # TODO: how many scenarios reduce to, i.e. sc
        #pro = scenario_reduced(pro, sc, scenario_reduction_clustering, 60)
        pro = create_stochastic_problem(problem_det_noise, sc, n_cargo_unknownweight[1], [], Bootstrap_bookedweight_quantile)
        # Save problem
        foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)_slack_$(slack_fraction[parse_index])"
        write_problem_stochastic(pro, foldername, "Stochastic_Problem", HPC_folder)
        #mo = create_model_stochastic(pro)
        mo = create_model_stochastic_slack_deck(pro, slack_fraction[parse_index]) 
        set_silent(mo) # removes terminal output
        set_time_limit_sec(mo, time_limit)
        optimize!(mo)
        sol = extract_stochastic_solution(pro, mo)
        # Save solution
        write_solution_stochastic(sol, foldername, "Stochastic_Solution", HPC_folder)
        cs_sol = sol.cs
        # Solve second stage when we know unknown weights
        #second_stage_m = second_stage_model(cs_sol, problem_det)
        second_stage_m = second_stage_model_deckslacked(cs_sol, problem_det, slack_fraction[parse_index])
        set_silent(second_stage_m) # removes terminal output
        set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
        optimize!(second_stage_m)
        fitted_sol = get_solution_second_stage_stochastic(problem_det, second_stage_m, sol)
        # Save fitted solution
        write_solution(fitted_sol, foldername, "Fitted_Solution", HPC_folder)
        # Check if second-stage problem was feasible
        #if fitted_sol.status == "INFEASIBLE"
        if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
            #second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
            second_stage_m_slacked = second_stage_model_deckslacked_slack(cs_sol, problem_det, slack_fraction[parse_index])
            set_silent(second_stage_m_slacked) # removes terminal output
            set_time_limit_sec(second_stage_m_slacked, time_limit)
            optimize!(second_stage_m_slacked)
            fitted_sol_slacked = get_solution_second_stage_stochastic(problem_det, second_stage_m_slacked, sol)
            write_solution(fitted_sol_slacked, foldername, "Fitted_Solution_slacked", HPC_folder)
            if primal_status(second_stage_m_slacked) == MOI.FEASIBLE_POINT || termination_status(second_stage_m_slacked) == MOI.OPTIMAL
                write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
            end
        end
        # EVP method for Bootstrap method 1
        foldername = "EVP_Bootstrap1_rep$(i)_sc$(sc)_unknown$(n_cargo_unknownweight[1])_time$(time_limit)_slack_$(slack_fraction[parse_index])"
        pro_sto = pro
        pro = expected_value_problem(pro)
        # Save problem
        write_problem(pro, foldername, "EVP_Problem", HPC_folder)
        #mo = create_model(pro)
        mo = create_model_slack_deck(pro, slack_fraction[parse_index])
        set_silent(mo) # removes terminal output
        set_time_limit_sec(mo, time_limit)
        optimize!(mo)
        sol = extract_solution(pro, mo)
        # Save solution
        write_solution(sol, foldername, "EVP_Solution", HPC_folder)
        cs_sol = sol.cs
        # Finding solutions to each scenario for EVP
        for l in 1:sc
            # create problems
            pro_temp = StowageProblem(
                vessel = problem_det.vessel,
                slots = problem_det.slots,
                cargo = pro_sto.cargo.items[l],    
                # Problem metadata
                name = "EVP Bootstrap $(l)",
                timestamp = now()
            )
            EVP_model = second_stage_model_deckslacked_slack(cs_sol, pro_temp, slack_fraction[parse_index])
            set_silent(EVP_model) # removes terminal output
            set_time_limit_sec(EVP_model, time_limit)
            optimize!(EVP_model)
            fitted_sol_EVP = get_solution_second_stage_deterministic(pro_temp, EVP_model, sol)
            #write_solution(fitted_sol_EVP, foldername, "Fitted_Solution_EVP_$(l)", HPC_folder)
            # Check if a feasible solution was found
            if primal_status(EVP_model) == MOI.FEASIBLE_POINT || termination_status(EVP_model) == MOI.OPTIMAL
                slack = [value.(EVP_model[:slack_deck]), value.(EVP_model[:slack_Vmax]),value.(EVP_model[:slack_Vmin]),
                    value.(EVP_model[:slack_Tmin]), value.(EVP_model[:slack_Tmax]), 
                    value.(EVP_model[:slack_Lmin]), value.(EVP_model[:slack_Lmax]),
                    #value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
                    value.(EVP_model[:slack_shearMin]),
                    value.(EVP_model[:slack_shearMax]), value.(EVP_model[:slack_bendingMax]),
                    value.(EVP_model[:slack_ballast_tanks])
                    ]
                # Check if slack is used
                if any(x->x > 0.001,vcat(slack...))
                    if sum(slack[1]) > 0.001 # deck limit
                        inf_index_EVP_boot[i,l] = 1
                    else # something else
                        inf_index_EVP_boot[i,l] = 2
                    end
                end
            else
                inf_index_EVP_boot[i,l] = -1
            end
            obj_val_EVP_boot[i,l] = fitted_sol_EVP.objective
            cargo_load_EVP_boot[i,l] = fitted_sol_EVP.n_cargo_loaded
            ballast_used_EVP_boot[i,l] = fitted_sol_EVP.ballast_weight   
        end

        # Second stage model with this solution and deterministic weights
        #second_stage_m = second_stage_model(cs_sol, problem_det)
        #set_silent(second_stage_m) # removes terminal output
        #set_time_limit_sec(second_stage_m, time_limit) # 5 minutes to solve model
        #optimize!(second_stage_m)
        #fitted_sol = get_solution_second_stage_deterministic(problem_det, second_stage_m, sol)
        # Save fitted solution
        #write_solution(fitted_sol, foldername, "Fitted_Solution", HPC_folder)
        # Check if second-stage problem was feasible
        #fitted_sol.status == "INFEASIBLE"
        #if fitted_sol.status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
            #second_stage_m_slacked = second_stage_model_slack(cs_sol, problem_det)
        #    second_stage_m_slacked = second_stage_model_deckslacked_slack(cs_sol, problem_det, slack_fraction[j])
        #    set_silent(second_stage_m_slacked) # removes terminal output
        #    set_time_limit_sec(second_stage_m_slacked, time_limit)
        #    optimize!(second_stage_m_slacked)
        #    fitted_sol_slacked = get_solution_second_stage_deterministic(problem_det, second_stage_m_slacked, sol)
        #    write_solution(fitted_sol_slacked, foldername, "Fitted_Solution_slacked", HPC_folder)
        #    if primal_status(second_stage_m_slacked) == MOI.FEASIBLE_POINT || termination_status(second_stage_m_slacked) == MOI.OPTIMAL
        #        write_slack(HPC_folder, foldername, "Fitted_Solution", second_stage_m_slacked)
        #    end
        #end
    end
#end

# Save results 
# gen
save_EVP_objective_2D(HPC_folder,"EVP_gen_info_slack_$(slack_fraction[parse_index])", "Objective_values_gen", obj_val_EVP_gen)
save_EVP_objective_2D(HPC_folder,"EVP_gen_info_slack_$(slack_fraction[parse_index])", "Cargo_loaded_gen", cargo_load_EVP_gen)
save_EVP_objective_2D(HPC_folder,"EVP_gen_info_slack_$(slack_fraction[parse_index])", "Ballast_used_gen", ballast_used_EVP_gen)
save_EVP_objective_2D(HPC_folder,"EVP_gen_info_slack_$(slack_fraction[parse_index])", "inf_index_gen", inf_index_EVP_gen)
# boot
save_EVP_objective_2D(HPC_folder,"EVP_boot_info_slack_$(slack_fraction[parse_index])", "Objective_values_boot", obj_val_EVP_boot)
save_EVP_objective_2D(HPC_folder,"EVP_boot_info_slack_$(slack_fraction[parse_index])", "Cargo_loaded_boot", cargo_load_EVP_boot)
save_EVP_objective_2D(HPC_folder,"EVP_boot_info_slack_$(slack_fraction[parse_index])", "Ballast_used_boot", ballast_used_EVP_boot)
save_EVP_objective_2D(HPC_folder,"EVP_boot_info_slack_$(slack_fraction[parse_index])", "inf_index_boot", inf_index_EVP_boot)
