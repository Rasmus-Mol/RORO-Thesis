# Main for stability test for stochastic solutions
include("packages_and_files.jl")

HPC_folder = "Finlandia_stability_tests_stochastic"
problemname1, problemname3 = "finlandia", "hazardous"
no_test = 10
repetitions = 10

# A problem for each instance for each repetition run
Sto_pro = Array{Any}(undef, 8, repetitions)

slack_sol_gen = Array{Any}(nothing,8, repetitions, no_test)
sol_gen = Array{Any}(undef, 8, repetitions, no_test)


for i in 1:8
    test_problem_name = Finlandia_test[i]
    for j in 1:repetitions
        # Load the problem
        foldername = "Stochastic_Stability_randomscenarios_uniformsampling_instance_$(i)_repetition_$(j)"
        filename = "Stochastic_Problem"
        Sto_pro[i,j] = get_stochastic_problem(foldername, filename, HPC_folder, 
                        problemname1, test_problem_name, problemname3)
        # Load the solution
        path = joinpath("Results", HPC_folder,foldername)
        files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
        for k in 1:no_test
            if "Solution_$(k)_info.json" in files # Could be solved
                sol_gen[i,j,k] = get_solution_deterministic(foldername,"Solution_$(k)",HPC_folder)
            else # was infeasible and therefore slacked
                sol_gen[i,j,k] = get_solution_deterministic(foldername,"Solution_$(k)_slacked",HPC_folder)
                # Should work when next test has been run 
                slack_sol_gen[i,j,k]  = get_slack(foldername,"Fitted_Solution_slacked_$(k)",HPC_folder)
                #slack_sol_gen_idx[i,j] = 1
            end
        end

    end
end

inf_count_gen = zeros(8,repetitions)
# count how many of the test were feasible
for i in 1:8 # number of instances
    for j in 1:repetitions # number of repetitions
        for k in 1:no_test
            if !isnothing(slack_sol_gen[i,j,k]) # if there is a slack solution
                inf_count_gen[i,j] += 1 # count as infeasible
                println("instance: $(i), Repetition: $(j), Test: $(k) was infeasible")
                if sum(slack_sol_gen[i,j,k][1]) < 1 # if there is a slack
                    println("Not deck 1 weight limit that is the issue: ", slack_sol_gen[i,j,k])
                end
            end
        end
    end
end
inf_count_gen

slack_sol_gen[8,1,1]
sol_gen[1,1,1].ballast_weight
sol_gen[1,1,4].ballast_weight