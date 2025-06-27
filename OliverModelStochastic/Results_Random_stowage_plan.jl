include("packages_and_files.jl")

problemname1, problemname3 = "finlandia", "hazardous"

# new results - old results are on USB
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
Det_sol = Array{Any}(undef, length(HPC_folders))
Det_pro = Array{Any}(undef, length(HPC_folders))
for l in 1:length(HPC_folders)
    test_instance = Finlandia_test[l]
    Det_sol[l] = get_solution_deterministic("Finlandia_deterministic",
        "Deterministic_Solution", HPC_folders[l])
    Det_pro[l] = load_data(problemname1, test_instance, problemname3)
    
end


HPC_folders = [
    "Finlandia_mixed_light_60_15_05_13",
    "Finlandia_mixed_light_100_15_05_15",
    "Finlandia_mixed_heavy_60_15_05_17",
    "Finlandia_mixed_heavy_100_15_05_17",
    "Finlandia_no_cars_light_60_14_05_20",
    "Finlandia_no_cars_light_100_15_05_10",
    "Finlandia_no_cars_heavy_60_15_05_10",
    "Finlandia_no_cars_heavy_100_15_05_09",
]

no_test1 = 10
no_test2 = 10
plot_folder = "Plots/Results/Random_plans/"

# Load problems and results
no_folders = length(HPC_folders)
test1_problems_gen = Array{Any}(undef, no_folders)
test1_problems_boot = Array{Any}(undef, no_folders)
test1_solutions_gen = Array{Any}(undef, no_folders,no_test1)
test1_solutions_boot = Array{Any}(undef, no_folders,no_test1)
slack_sol_gen = Array{Any}(undef, no_folders,no_test1)
slack_sol_boot = Array{Any}(undef, no_folders,no_test1)
slack_sol_gen_idx = zeros(no_folders,no_test1)
slack_sol_boot_idx = zeros(no_folders,no_test1)

initial_random_plan_trucksfirst = Array{Any}(undef, no_folders,no_test2)
initial_random_plan_carsfirst = Array{Any}(undef, no_folders,no_test2)
initial_random_plan_secufirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_trucksfirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_carsfirst = Array{Any}(undef, no_folders,no_test2)
sol_no_shifts_secufirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_trucksfirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_carsfirst = Array{Any}(undef, no_folders,no_test2)
sol_shifts_secufirst = Array{Any}(undef, no_folders,no_test2)

Det_sol = Array{Any}(undef, no_folders)
Det_pro = Array{Any}(undef, no_folders)

for i in 1:no_folders
    test_problem_name = Finlandia_test[i]
    HPC_folder_load = HPC_folders[i]
    # load deterministic solution and problem
    Det_sol[i] =  get_solution_deterministic("Finlandia_deterministic",
    "Deterministic_Solution", HPC_folder_load)
    Det_pro[i] = load_data(problemname1, test_problem_name, problemname3)
    # Test 1 - Uniform sampling
    foldername = "Determinitic_Stability_randomscenarios_uniformsampling"
    filename = "Stochastic_Problem"
    test1_problems_gen[i] = get_stochastic_problem(foldername,
                filename, HPC_folder_load, problemname1, test_problem_name, problemname3)
    path = joinpath("Results", HPC_folder_load,foldername)
    files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
    for j in 1:no_test1
        if "Solution_$(j)_info.json" in files # Could be solved
            test1_solutions_gen[i,j] = get_solution_deterministic(foldername,"Solution_$(j)",HPC_folder_load)
        else # was infeasible and therefore slacked
            test1_solutions_gen[i,j] = get_solution_deterministic(foldername,"Solution_$(j)_slacked",HPC_folder_load)
            # Should work when next test has been run 
            slack_sol_gen[i,j]  = get_slack(foldername,"Fitted_Solution_slacked_$(j)",HPC_folder_load)
            slack_sol_gen_idx[i,j] = 1
        end
    end
    foldername = "Determinitic_Stability_randomscenarios_bootstrapsampling"
    filename = "Stochastic_Problem"
    test1_problems_boot[i] = get_stochastic_problem(foldername,
                filename, HPC_folder_load, problemname1, test_problem_name, problemname3)
    path = joinpath("Results", HPC_folder_load,foldername)
    files = filter(f -> startswith(f, "Solution_") && isfile(joinpath(path, f)), readdir(path))
    for j in 1:no_test1
        if "Solution_$(j)_info.json" in files # Could be solved
            test1_solutions_boot[i,j] = get_solution_deterministic(foldername,"Solution_$(j)",HPC_folder_load)
        else # was infeasible and therefore slacked
            test1_solutions_boot[i,j] = get_solution_deterministic(foldername,"Solution_$(j)_slacked",HPC_folder_load)
            # Should work when next test has been run 
            slack_sol_boot[i,j] = get_slack(foldername,"Fitted_Solution_slacked_$(j)",HPC_folder_load)
            slack_sol_boot_idx[i,j] = 1
        end
    end
    # Test 2 - random plan
    for j in 1:no_test2
        # Trucks first
        foldername = "Random_Stowage_Plan_Trucks_first"
        #filename = "Random_plan" # should be fixed next time
        filename = "Random_plan_$(j)" 
        initial_random_plan_trucksfirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_trucksfirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(j)",HPC_folder_load)
        sol_shifts_trucksfirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(j)",HPC_folder_load)
        # Cars first
        foldername = "Random_Stowage_Plan_Cars_first"
        #filename = "Random_plan" # should be fixed next time
        filename = "Random_plan_$(j)" 
        initial_random_plan_carsfirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_carsfirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(j)",HPC_folder_load)
        sol_shifts_carsfirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(j)",HPC_folder_load)
        # Secu first
        foldername = "Random_Stowage_Plan_Secu_first"
        #filename = "Random_plan" # should be fixed next time
        filename = "Random_plan_$(j)" 
        initial_random_plan_secufirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
        sol_no_shifts_secufirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(j)",HPC_folder_load)
        sol_shifts_secufirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(j)",HPC_folder_load) 
    end

end
# Empty ship
HPC_folder_load = HPC_folders[1]
sol_empty_ship = get_solution_deterministic("Finlandia_deterministic","Empty_ship", HPC_folder_load)
# Empty ship and no ballast water
empty_ship_no_ballast = empty_ship_model(Det_pro[1])
set_silent!(empty_ship_no_ballast)
optimize!(empty_ship_no_ballast)
# check if model is infeasible
if termination_status(empty_ship_no_ballast) in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
    println("Empty ship model is infeasible")
else
    println("Empty ship model is feasible")
end
# Why is ship infeasible


#######################
# test 1
idx_gen = findall(x -> x==1, slack_sol_gen_idx)
idx_boot = findall(x -> x==1, slack_sol_boot_idx)
println("Test that were infeasible - gen: ",idx_gen)
println("Test that were infeasible - boot: ",idx_boot)
println("number of test that were infeasible: ", Int(sum(slack_sol_gen_idx)),"/", no_folders*no_test1)
println("number of test that were infeasible boot: ", Int(sum(slack_sol_boot_idx)),"/", no_folders*no_test1)
# Problems that were still feasible
for i in 1:no_folders
    println("Deterministic total weight: ",Det_pro[i].cargo.total_weight)
    for j in 1:no_test1
        if slack_sol_gen_idx[i,j] == 0 # was feasible
            println("scenario total weight: ",test1_problems_gen[i].cargo.items[j].total_weight)
        end
    end
end
# Change in ballast water
ballast_water_gen = zeros(no_folders,no_test1)
ballast_water_boot = zeros(no_folders,no_test1)
avg_ballast_water_gen = zeros(no_folders)
avg_ballast_water_boot = zeros(no_folders)
for i in 1:no_folders
    for j in 1:no_test1
        if slack_sol_gen_idx[i,j] == 0 # was feasible
            ballast_water_gen[i,j] = test1_solutions_gen[i,j].ballast_weight
        end
        if slack_sol_boot_idx[i,j] == 0 # was feasible
            ballast_water_boot[i,j] = test1_solutions_boot[i,j].ballast_weight
        end
    end
    println("####################")
    println("Deterministic ballast water: ", Det_sol[i].ballast_weight)
    println("Ballast water gen: ", ballast_water_gen[i,:])
    println("Ballast water boot: ", ballast_water_boot[i,:])
    avg_ballast_water_boot[i] = mean(filter(x-> x>0,ballast_water_boot[i,:]))
    avg_ballast_water_gen[i] = mean(filter(x-> x>0,ballast_water_gen[i,:]))
end
for i in 1:no_folders
    #println(round(Det_sol[i].ballast_weight,digits = 2))
    #println(round(Det_pro[i].cargo.total_weight,digits = 2))
    #println(sum(Det_pro[i].cargo.items[j].))
    #println(round(mean(test1_problems_boot[i].cargo.items.total_weight),digits = 2))
    #println(round(mean(test1_problems_gen[i].cargo.items.total_weight),digits = 2))
    #println(no_test1-sum(slack_sol_gen_idx[i,:]))
    #println(no_test1-sum(slack_sol_boot_idx[i,:]))
    #println(round(avg_ballast_water_gen[i],digits=2))
    println(round(avg_ballast_water_boot[i],digits=2))
end

# Slack order:
# slack_deck, slack_Vmax, slack_Vmin, slack_Tmin, slack_Tmax,
# slack_Lmin, slack_Lmax, slack_shearMin, slack_shearMax, slack_bendingMax,
# slack_ballast_tanks
# infeasible because: Needs new results before doing
inf_gen = []
inf_boot = []
id_deck1_gen = Matrix{Vector{Int64}}(undef,no_folders,no_test1)
id_deck1_boot = Matrix{Vector{Int64}}(undef,no_folders,no_test1)
weight_deck_1_gen = Matrix{Float64}(undef,no_folders,no_test1)
weight_deck_1_boot = Matrix{Float64}(undef,no_folders,no_test1)

id_deck2_gen = Matrix{Vector{Int64}}(undef,no_folders,no_test1)
id_deck2_boot = Matrix{Vector{Int64}}(undef,no_folders,no_test1)
weight_deck_2_gen = Matrix{Float64}(undef,no_folders,no_test1)
weight_deck_2_boot = Matrix{Float64}(undef,no_folders,no_test1)

id_deck3_gen = Matrix{Vector{Int64}}(undef,no_folders,no_test1)
id_deck3_boot = Matrix{Vector{Int64}}(undef,no_folders,no_test1)
weight_deck_3_gen = Matrix{Float64}(undef,no_folders,no_test1)
weight_deck_3_boot = Matrix{Float64}(undef,no_folders,no_test1)

for i in 1:no_folders
    for j in 1:no_test1
        id_deck1_gen[i,j] = [c.id for c in filter(x -> x.deck==1,test1_solutions_gen[i,j].cargo)]
        id_deck1_boot[i,j] = [c.id for c in filter(x -> x.deck==1,test1_solutions_boot[i,j].cargo)]
        weight_deck_1_gen[i,j] = sum([c.weight for c in test1_problems_gen[i].cargo.items[j] if c.id in id_deck1_gen[i,j]])
        weight_deck_1_boot[i,j] = sum([c.weight for c in test1_problems_boot[i].cargo.items[j] if c.id in id_deck1_boot[i,j]])
        #weight_deck_1_gen[i,j] = sum([c.weight for c in filter(x -> x.deck==1,test1_solutions_gen[i,j].cargo)])
        #weight_deck_1_boot[i,j] = sum([c.weight for c in filter(x -> x.deck==1,test1_solutions_boot[i,j].cargo)])
        id_deck2_gen[i,j] = [c.id for c in filter(x -> x.deck==2,test1_solutions_gen[i,j].cargo)]
        id_deck2_boot[i,j] = [c.id for c in filter(x -> x.deck==2,test1_solutions_boot[i,j].cargo)]
        weight_deck_2_gen[i,j] = sum([c.weight for c in test1_problems_gen[i].cargo.items[j] if c.id in id_deck2_gen[i,j]])
        weight_deck_2_boot[i,j] = sum([c.weight for c in test1_problems_boot[i].cargo.items[j] if c.id in id_deck2_boot[i,j]])
        id_deck3_gen[i,j] = [c.id for c in filter(x -> x.deck==3,test1_solutions_gen[i,j].cargo)]
        id_deck3_boot[i,j] = [c.id for c in filter(x -> x.deck==3,test1_solutions_boot[i,j].cargo)]
        weight_deck_3_gen[i,j] = sum([c.weight for c in test1_problems_gen[i].cargo.items[j] if c.id in id_deck3_gen[i,j]])
        weight_deck_3_boot[i,j] = sum([c.weight for c in test1_problems_boot[i].cargo.items[j] if c.id in id_deck3_boot[i,j]])
        if slack_sol_gen_idx[i,j] == 1 # was infeasible
            push!(inf_gen,findall(x-> x>0,[sum(slack_sol_gen[i,j][k]) for k in 1:length(slack_sol_gen[i,j])]))
            if inf_gen[end][1] == 1
                println("Deck limit: ", slack_sol_gen[i,j][1])
            end
        end
        if slack_sol_boot_idx[i,j] == 1 # was infeasible
            push!(inf_boot,findall(x-> x>0,[sum(slack_sol_boot[i,j][k]) for k in 1:length(slack_sol_boot[i,j])]))
        end
    end
end
# deck 1 plots
p1 = plot(weight_deck_1_gen[1,:],label = "Instance 1",xlabel = "Test", ylabel = "Weight (t.)",
    title = "Weight on deck 1 - SGM. 1")
for i in 2:no_folders
    plot!(p1,weight_deck_1_gen[i,:],label = "Instance $(i)")
end
plot!(ones(no_test1)*Det_pro[1].vessel.decks[1].weight_limit, label = "Deck 1 limit", color = :red, linestyle = :dash)
display(p1)
savefig(p1, plot_folder*"Weight_Deck1_SGM1.png")
p2 = plot(weight_deck_1_boot[1,:],label = "Instance 1",xlabel = "Test", ylabel = "Weight (t.)",
    title = "Weight on deck 1 - SGM. 2")
for i in 2:no_folders
    plot!(p2,weight_deck_1_boot[i,:],label = "Instance $(i)")
end
plot!(p2,ones(no_test1)*Det_pro[1].vessel.decks[1].weight_limit, label = "Deck 1 limit", color = :red, linestyle = :dash)
savefig(p2, plot_folder*"Weight_Deck1_SGM2.png")
display(p2)
det_sol_weight_deck1 = []
no_folders = length(HPC_folders)
for i in 1:no_folders
    push!(det_sol_weight_deck1, sum([c.weight for c in filter(x -> x.deck==1,Det_sol[i].cargo)]))
end
p3 = plot(det_sol_weight_deck1,xlabel = "Instance", ylabel = "Weight (t.)",
    title = "Weight on deck 1 - Deterministic solution", label = "Deterministic solution")
plot!(p3, ones(8)*Det_pro[1].vessel.decks[1].weight_limit, label = "Deck 1 limit", color = :red, linestyle = :dash)
savefig(p3, plot_folder*"Weight_Deck1_Deterministic.png")
# deck 2 plots
p1 = plot(weight_deck_2_gen[1,:],label = "Instance 1",xlabel = "Test", ylabel = "Weight (t.)",
    title = "Weight on deck 2 - SGM. 1")
for i in 2:no_folders
    plot!(p1,weight_deck_2_gen[i,:],label = "Instance $(i)")
end
plot!(ones(no_test1)*Det_pro[1].vessel.decks[2].weight_limit, label = "Deck 2 limit", color = :red, linestyle = :dash)
display(p1)
savefig(p1, plot_folder*"Weight_Deck2_SGM1.png")
p2 = plot(weight_deck_2_boot[1,:],label = "Instance 1",xlabel = "Test", ylabel = "Weight (t.)",
    title = "Weight on deck 2 - SGM. 2")
for i in 2:no_folders
    plot!(p2,weight_deck_2_boot[i,:],label = "Instance $(i)")
end
plot!(p2,ones(no_test1)*Det_pro[1].vessel.decks[2].weight_limit, label = "Deck 2 limit", color = :red, linestyle = :dash)
savefig(p2, plot_folder*"Weight_Deck2_SGM2.png")
display(p2)
det_sol_weight_deck2 = []
for i in 1:no_folders
    push!(det_sol_weight_deck2, sum([c.weight for c in filter(x -> x.deck==2,Det_sol[i].cargo)]))
end
p3 = plot(det_sol_weight_deck2,xlabel = "Instance", ylabel = "Weight (t.)",
    title = "Weight on deck 2 - Deterministic solution", label = "Deterministic solution")
plot!(p3, ones(8)*Det_pro[1].vessel.decks[2].weight_limit, label = "Deck 2 limit", color = :red, linestyle = :dash)
savefig(p3, plot_folder*"Weight_Deck2_Deterministic.png")
display(p3)
# deck 3 plots
p1 = plot(weight_deck_3_gen[1,:],label = "Instance 1",xlabel = "Test", ylabel = "Weight (t.)",
    title = "Weight on deck 3 - SGM. 1")
for i in 2:no_folders
    plot!(p1,weight_deck_3_gen[i,:],label = "Instance $(i)")
end
plot!(ones(8)*Det_pro[1].vessel.decks[3].weight_limit, label = "Deck 3 limit", color = :red, linestyle = :dash)
display(p1)
savefig(p1, plot_folder*"Weight_Deck3_SGM1.png")
p2 = plot(weight_deck_3_boot[1,:],label = "Instance 1",xlabel = "Test", ylabel = "Weight (t.)",
    title = "Weight on deck 3 - SGM. 2")
for i in 2:no_folders
    plot!(p2,weight_deck_3_boot[i,:],label = "Instance $(i)")
end
plot!(p2,ones(8)*Det_pro[1].vessel.decks[3].weight_limit, label = "Deck 3 limit", color = :red, linestyle = :dash)
savefig(p2, plot_folder*"Weight_Deck3_SGM2.png")
display(p2)
det_sol_weight_deck3 = []
for i in 1:no_folders
    push!(det_sol_weight_deck3, sum([c.weight for c in filter(x -> x.deck==3,Det_sol[i].cargo)]))
end
p3 = plot(det_sol_weight_deck3,xlabel = "Instance", ylabel = "Weight (t.)",
    title = "Weight on deck 3 - Deterministic solution", label = "Deterministic solution")
plot!(p3, ones(8)*Det_pro[1].vessel.decks[3].weight_limit, label = "Deck 3 limit", color = :red, linestyle = :dash)
savefig(p3, plot_folder*"Weight_Deck3_Deterministic.png")
det_sol_weight_deck1[4]
Det_pro[1].vessel.decks[1].weight_limit
weight_cars = []
weight_trucks = []
weight_secus = []
weight_heavy_machinery = []
weight_cars_gen = Matrix{Float64}(undef,no_folders,no_test1)
weight_trucks_gen = Matrix{Float64}(undef,no_folders,no_test1)
weight_secus_gen = Matrix{Float64}(undef,no_folders,no_test1)
weight_heavy_machinery_gen = Matrix{Float64}(undef,no_folders,no_test1)
weight_cars_boot = Matrix{Float64}(undef,no_folders,no_test1)
weight_trucks_boot = Matrix{Float64}(undef,no_folders,no_test1)
weight_secus_boot = Matrix{Float64}(undef,no_folders,no_test1)
weight_heavy_machinery_boot = Matrix{Float64}(undef,no_folders,no_test1)
# Collect weight of each cargo type
for i in 1:4
    # Det
    push!(weight_cars, sum([c.weight for c in filter(x -> x.cargo_type_id == 2,Det_pro[i].cargo)]))
    push!(weight_trucks, sum([c.weight for c in filter(x -> x.cargo_type_id == 1,Det_pro[i].cargo)]))
    push!(weight_secus, sum([c.weight for c in filter(x -> x.cargo_type_id == 4,Det_pro[i].cargo)]))
    push!(weight_heavy_machinery, sum([c.weight for c in filter(x -> x.cargo_type_id == 3,Det_pro[i].cargo)]))
    for j in 1:no_test1
        # Gen
        weight_cars_gen[i,j]= sum([c.weight for c in filter(x -> x.cargo_type_id == 2,test1_problems_gen[i].cargo.items[j])])
        weight_trucks_gen[i,j]= sum([c.weight for c in filter(x -> x.cargo_type_id == 1,test1_problems_gen[i].cargo.items[j])])
        weight_secus_gen[i,j]= sum([c.weight for c in filter(x -> x.cargo_type_id == 4,test1_problems_gen[i].cargo.items[j])])
        weight_heavy_machinery_gen[i,j]= sum([c.weight for c in filter(x -> x.cargo_type_id == 3,test1_problems_gen[i].cargo.items[j])])
        # Boot
        weight_cars_boot[i,j]= sum([c.weight for c in filter(x -> x.cargo_type_id == 2,test1_problems_boot[i].cargo.items[j])])
        weight_trucks_boot[i,j]= sum([c.weight for c in filter(x -> x.cargo_type_id == 1,test1_problems_boot[i].cargo.items[j])])
        weight_secus_boot[i,j]= sum([c.weight for c in filter(x -> x.cargo_type_id == 4,test1_problems_boot[i].cargo.items[j])])
        weight_heavy_machinery_boot[i,j]= sum([c.weight for c in filter(x -> x.cargo_type_id == 3,test1_problems_boot[i].cargo.items[j])])
    end
end
# cars
println(weight_cars[4])
println(weight_cars_gen[4,:])
# Trucks
println(weight_trucks[4])
println(weight_trucks_gen[4,:])
# Secus
println(weight_secus[4])
println(weight_secus_gen[4,:])
# Heavy machinery
println(weight_heavy_machinery[4])
println(weight_heavy_machinery_gen[4,:])

######
pl1 = Det_sol[1].cargo
sum([c.weight for c in filter(x -> x.deck==1,pl1)])
length([c.weight for c in filter(x -> x.deck==1,pl1)])
length([c.weight for c in filter(x -> x.deck==2,pl1)])
length([c.weight for c in filter(x -> x.deck==3,pl1)])
pl2 = Det_sol[2].cargo
sum([c.weight for c in filter(x -> x.deck==1,pl2)])
length([c.weight for c in filter(x -> x.deck==1,pl2)])
length([c.weight for c in filter(x -> x.deck==2,pl2)])
length([c.weight for c in filter(x -> x.deck==3,pl2)])


p1 = plot_solution(Det_sol[1])
p2 = plot_solution(Det_sol[2])
p3 = plot_solution(Det_sol[3])
p4 = plot_solution(Det_sol[4])
p5 = plot_solution(Det_sol[5])
p6 = plot_solution(Det_sol[6])
p7 = plot_solution(Det_sol[7])
p8 = plot_solution(Det_sol[8])
savefig(p1, plot_folder*"Deterministic_solution_instance_1.png")
savefig(p2, plot_folder*"Deterministic_solution_instance_2.png")
savefig(p3, plot_folder*"Deterministic_solution_instance_3.png")
savefig(p4, plot_folder*"Deterministic_solution_instance_4.png")
savefig(p5, plot_folder*"Deterministic_solution_instance_5.png")
savefig(p6, plot_folder*"Deterministic_solution_instance_6.png")
savefig(p7, plot_folder*"Deterministic_solution_instance_7.png")
savefig(p8, plot_folder*"Deterministic_solution_instance_8.png")


#################################

# test 2 - Stability of random plan
shifts_carsfirst = zeros(no_folders,no_test2)
shifts_trucksfirst = zeros(no_folders,no_test2)
shifts_secufirst = zeros(no_folders,no_test2)
random_plan_sol_truck = Array{Any}(nothing, no_folders,no_test2)
random_plan_sol_cars = Array{Any}(nothing, no_folders,no_test2)
random_plan_sol_secu = Array{Any}(nothing, no_folders,no_test2)
slack_truck = Array{Any}(nothing, no_folders,no_test2)
slack_cars = Array{Any}(nothing, no_folders,no_test2)
slack_secu = Array{Any}(nothing, no_folders,no_test2)
for i in 1:no_folders
    for j in 1:no_test2
        shifts_trucksfirst[i,j] = count((!=).(initial_random_plan_trucksfirst[i,j].cs,
                                sol_no_shifts_trucksfirst[i,j].cs))
        random_plan_sol_truck[i,j] = get_solution_random_stowage_plan(initial_random_plan_trucksfirst[i,j],Det_pro[i])
        if shifts_trucksfirst[i,j] >0
            slack_truck[i,j] = model_with_slack_get_slac(initial_random_plan_trucksfirst[i,j],Det_pro[i])
        end

        shifts_carsfirst[i,j] = count((!=).(initial_random_plan_carsfirst[i,j].cs,
                                sol_no_shifts_carsfirst[i,j].cs))
        random_plan_sol_cars[i,j] = get_solution_random_stowage_plan(initial_random_plan_carsfirst[i,j],Det_pro[i])
        if shifts_carsfirst[i,j] >0
            slack_cars[i,j] = model_with_slack_get_slac(initial_random_plan_carsfirst[i,j],Det_pro[i])
        end

        shifts_secufirst[i,j] = count((!=).(initial_random_plan_secufirst[i,j].cs,
                                sol_no_shifts_secufirst[i,j].cs))
        random_plan_sol_secu[i,j] = get_solution_random_stowage_plan(initial_random_plan_secufirst[i,j],Det_pro[i])
        if shifts_secufirst[i,j] >0
            slack_secu[i,j] = model_with_slack_get_slac(initial_random_plan_secufirst[i,j],Det_pro[i])  
        end
    end 
end

for i in 1:no_folders
    #println("instance $(j): ",mean([sum(initial_random_plan_trucksfirst[j,i].cs) for i in 1:no_test2]))
    #println("instance $(j): ",mean([sum(initial_random_plan_carsfirst[j,i].cs) for i in 1:no_test2]))
    #println("instance $(j): ",mean([sum(initial_random_plan_secufirst[j,i].cs) for i in 1:no_test2]))
    #id_shift_trucks = findall(x -> x>0, shifts_trucksfirst[i,:])
    #id_shift_cars = findall(x -> x>0, shifts_carsfirst[i,:])
    #id_shift_secu = findall(x -> x>0, shifts_secufirst[i,:])
    #println("Number of times random plan was infeasible, instance $(i):")
    #println("Trucks First: ", length(id_shift_trucks))
    #println("Cars First: ", length(id_shift_cars))
    #println("Secu-boxes First: ", length(id_shift_secu))
    #for j in 1:no_test2
    #    println(initial_random_plan_secufirst[i,j].cs == sol_no_shifts_secufirst[i,j].cs)
    #end
    
    temp = 0
    a = 0
    for j in 1:no_test2
        if isnothing(slack_truck[i,j])
            #isnothing(slack_cars[i,j])
            #isnothing(slack_secu[i,j])
            #
            #temp +=sol_no_shifts_trucksfirst[i,j].ballast_weight
            #sol_no_shifts_carsfirst[i,j].ballast_weight
            #sol_no_shifts_secufirst[i,j].ballast_weight

            #a+=1
        end
        temp += sol_shifts_trucksfirst[i,j].ballast_weight
        #sol_shifts_secufirst[i,j].ballast_weight
        a+= 1
    end
    if a>0
        println("Mean ballast weight, instance $(i): ", round(temp/a,digits=2))
    end
    #println("a: ",a)
    
    #=
    for j in 1:no_test2
        if isnothing(slack_secu[i,j])
            #println("Instance $(i), test $(j): No slack for secu boxes")
        else
            println("Instance $(i), test $(j): Slack for secu boxes: ", slack_secu[i,j])
        end
        if isnothing(slack_cars[i,j])
            #println("Instance $(i), test $(j): No slack for cars")
        else
            println("Instance $(i), test $(j): Slack for cars: ", slack_cars[i,j])
        end
        if isnothing(slack_truck[i,j])
            #println("Instance $(i), test $(j): No slack for trucks")
        else
            println("Instance $(i), test $(j): Slack for trucks: ", slack_truck[i,j])
        end
    end
    =#
end
diff_placement = []
for i in 1:no_test2
    println("################")
    #println("Shifts:",    sol_shifts_trucksfirst[1,i].ballast_weight)
    #println("No shifts: ",    sol_no_shifts_trucksfirst[1,i].ballast_weight)
    push!(diff_placement,count((!=).(sol_shifts_trucksfirst[1,i].cs,
    sol_no_shifts_trucksfirst[1,i].cs)))
    #println(sol_shifts_trucksfirst[1,i].cs == sol_no_shifts_trucksfirst[1,i].cs)
end
diff_placement
plot_solution_random_plan(sol_shifts_trucksfirst[1,2])
plot_solution_random_plan(sol_no_shifts_trucksfirst[1,2])
# testing if it is actually correct
cargoc = Det_pro[1].cargo
vessel = Det_pro[1].vessel
slots = Det_pro[1].slots
ballast_weight_test = []
for j in 1:no_test2
    stowed_cargo = Cargo[]
    count = 1
    for i in 1:length(cargoc)
        if sum(initial_random_plan_carsfirst[1,j].cs[i,:])>0
            #sum(initial_random_plan_trucksfirst[1,j].cs[i,:])>0
            # add cargo to cargoplacement
            #push!(stowed_cargo, cargoc[i])
            push!(stowed_cargo, Cargo(
                id = count,
                cargo_type_id = cargoc[i].cargo_type_id,
                weight = cargoc[i].weight,
                loading_port = cargoc[i].loading_port,
                discharge_port = cargoc[i].discharge_port,
                priority = cargoc[i].priority,
                requires_lashing = cargoc[i].requires_lashing,
                requires_ventilation = cargoc[i].requires_ventilation,
                hazardous = cargoc[i].hazardous,
                refers = cargoc[i].refers)
            )
            count+=1
        end
    end
    new_cargoc = CargoCollection(stowed_cargo)
    # Creates Deterministic problem and model
    problem_det = StowageProblem(
        vessel = vessel,
        slots = slots,
        cargo = new_cargoc,
        name = Det_pro[1].name,
        timestamp = now()
    )
    model_det = create_model(problem_det)
    set_silent(model_det) # removes terminal output
    set_time_limit_sec(model_det, 60 * 5) # 1 hour
    optimize!(model_det)
    solution_det = extract_solution(problem_det, model_det)
    push!(ballast_weight_test, solution_det.ballast_weight)
end
println(sum(ballast_weight_test .- getfield.(sol_no_shifts_carsfirst[1,:],:ballast_weight)))
println(sum(ballast_weight_test .- getfield.(sol_shifts_carsfirst[1,:],:ballast_weight)))
#println(sum(ballast_weight_test .- getfield.(sol_no_shifts_trucksfirst[1,:],:ballast_weight)))
#println(sum(ballast_weight_test .- getfield.(sol_shifts_trucksfirst[1,:],:ballast_weight)))

p1 = plot_solution_random_plan(sol_shifts_trucksfirst[1,1])
savefig(p1, plot_folder*"Random_plan_solution__shifts_trucksfirst_1.png")
p2 = plot_solution_random_plan(sol_no_shifts_trucksfirst[1,1])
savefig(p2, plot_folder*"Random_plan_solution_no_shifts_trucksfirst_1.png")
#plot_solution(solution_det)
p3 = scatter(sol_shifts_trucksfirst[1,1].ballast_volume, label = "Shifts allowed. Sum: $(round(sum(sol_shifts_trucksfirst[1,1].ballast_volume),digits=2))", xlabel = "Ballast tank",
    ylabel = "Ballast weight (t.)", title = "Ballast weight in ballast tanks")
scatter!(p3, sol_no_shifts_trucksfirst[1,1].ballast_volume, label = "No shifts allowed. Sum: $(round(sum(sol_no_shifts_trucksfirst[1,1].ballast_volume),digits=2))", marker = :utriangle)
savefig(p3, plot_folder*"Ballast_weight_in_ballast_tanks.png")

slack_vars = Array{Any}(undef,no_test2)
slack_vars_string = Array{Any}(undef,no_test2)
for i in 1:no_test2
    s1 = []
    s2 = []
    #model_temp = No_ballast_slacked(sol_shifts_trucksfirst[1, i].cs, Det_pro[1])
    model_temp = No_ballast_slacked(sol_shifts_carsfirst[2, i].cs, Det_pro[2])
    set_silent(model_temp) # removes terminal output
    set_time_limit_sec(model_temp, 60 * 5) # 1 hour
    optimize!(model_temp)
    sol_temp = extract_solution(Det_pro[1], model_temp)
    slack = [value.(model_temp[:slack_deck]), value.(model_temp[:slack_Vmax]), value.(model_temp[:slack_Vmin]),
        value.(model_temp[:slack_Tmin]), value.(model_temp[:slack_Tmax]),
        value.(model_temp[:slack_Lmin]), value.(model_temp[:slack_Lmax]),
        #value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
        value.(model_temp[:slack_shearMin]),
        value.(model_temp[:slack_shearMax]), value.(model_temp[:slack_bendingMax]),
        value.(model_temp[:slack_ballast_tanks]),
        # New slack variable
        value.(model_temp[:slack_dis_min]), value.(model_temp[:slack_dis_max]),
        value.(model_temp[:slack_bending_pos]), value.(model_temp[:slack_bending_neg])
    ]
    tempmin = value.(collect(model_temp[:z_min]))
    tempmax = value.(collect(model_temp[:z_max]))
    println(typeof(tempmin))
    println("Z_min: ",findall(x->x>0,tempmin))
    println("Z_max: ",findall(x->x>0,tempmax))
    # string over slack variables:
    slack_string = ["slack deck", "slack Vmax", "slack Vmin", "slack Tmin", "slack Tmax",
        "slack Lmin", "slack Lmax", "slack shearMin", "slack shearMax",
        "slack bendingMax", "slack ballast tanks", "slack dis min",
        "slack dis max", "slack bending pos", "slack bending neg"]
    for j in 1:length(slack)
        if sum(slack[j]) > 0
            push!(s1, slack[j])
            push!(s2, slack_string[j])
        end
    end
    slack_vars[i] = s1
    slack_vars_string[i] = s2
end
for i in 1:no_test2
    println(slack_vars_string[i])
    for j in 1:length(slack_vars[i])
        println("Slack variable $(j): ", slack_vars[i][j])
    end
end
# Ship has to have a min sum of cargo and ballast weight
sol_empty_ship.ballast_weight+sol_empty_ship.cargo_weight
for i in 1:no_test2
    println("##############")
    println("test: ",i)
    println(round(sol_no_shifts_secufirst[1,i].cargo_weight,digits=2))
    println(round(sol_no_shifts_secufirst[1,i].ballast_weight,digits=2))
    println(round(sol_no_shifts_secufirst[1,i].cargo_weight+sol_no_shifts_secufirst[1,i].ballast_weight,digits=2))
    #println(sol_no_shifts_carsfirst[1,i].ballast_weight+sol_no_shifts_carsfirst[1,i].cargo_weight)
    #println(sol_no_shifts_carsfirst[2,i].ballast_weight+sol_no_shifts_carsfirst[2,i].cargo_weight)
end

displacement = [w.displacement for w in vessel.weight_dependencies]
sol_empty_ship.cargo_weight
# The cars first are interesting
sol_shifts_carsfirst[5,3].ballast_weight
sol_no_shifts_carsfirst[5,3].ballast_weight

plot_solution_random_plan(random_plan_sol_secu[8,2])
plot_solution_random_plan(sol_no_shifts_secufirst[8,2])
plot_solution_random_plan(random_plan_sol_truck[8,7])
plot_solution_random_plan(sol_no_shifts_trucksfirst[8,7])
plot_solution_random_plan(random_plan_sol_cars[8,7])
plot_solution_random_plan(sol_no_shifts_carsfirst[8,7])
n_cargo = length(Det_pro[8].cargo)
idx_not_stowed = []
for i in 1:n_cargo
    #if (sum(initial_random_plan_trucksfirst[4,1].cs[i,:]) < 0.5)
    #    push!(idx_not_stowed,i)
    #end
    if (sum(initial_random_plan_trucksfirst[8,1].cs[i,:]) < 0.5)
        push!(idx_not_stowed,i)
    end
end
idx_not_stowed


# Vi har random plan cs
# no shifts solution
# shift allowed solution
initial_random_plan_trucksfirst[i,j] = get_solution_deterministic(foldername,filename,HPC_folder_load)
sol_no_shifts_trucksfirst[i,j] = get_solution_deterministic(foldername,"Solution_no_shifts_$(j)",HPC_folder_load)
sol_shifts_trucksfirst[i,j] = get_solution_deterministic(foldername,"Solution_shifts_$(j)",HPC_folder_load)
        


# Takes random stowage plan, which only has the placement,
# and a problem, which also have weight of cargo, 
# and returns a solution with placement and weight
# so it can be plotted
function get_solution_random_stowage_plan(sol,pro)
    cs_val = sol.cs
    # Get placed cargo
    placements = CargoPlacement[] # Rasmus: Same as Vector{CargoPlacement}()
    n_cargo_loaded = 0
    for c in 1:length(pro.cargo), s in 1:length(pro.slots)
        if cs_val[c,s] > 0.5  # Binary variable threshold
            n_cargo_loaded += 1
            slot = pro.slots[s]
            cargo = pro.cargo[c]
            push!(placements, CargoPlacement(
                id = cargo.id,
                cargo_type_id = cargo.cargo_type_id,
                deck = slot.deck_id,
                lcg = slot.lcg,
                tcg = slot.tcg,
                vcg = slot.vcg,
                weight = cargo.weight,
                length = get_length(cargo),
                width = get_width(cargo),
                height = get_height(cargo),
                haz_class = cargo.hazardous
            ))
        end
    end
    return Solution(
        gap = Inf,
        status = 0,
        objective = Inf,
        time = 0,
        cargo_weight = 0.0,
        total_weight = 0.0,
        ballast_weight = 0.0, 
        area_utilization = 0.0,
        cargo = placements,
        cs = cs_val,
        n_cargo_total = sol.n_cargo_total,
        n_cargo_loaded = sol.n_cargo_loaded,
        shear_force = Float64[],
        bending_moment = Float64[],
        ballast_volume = Float64[],
        lcg = 0.0,
        tcg = 0.0, 
        vcg = 0.0,
        n_variables = 0,
        n_binary_variables = 0,
        n_constraints = 0,
        model_size = 0,
        solver_name = "Gurobi",
        solver_iterations = 0,
        solver_nodes = 0
    )
end

function model_random_plan_slacked(cs1, problem::StowageProblem)
    @unpack vessel, slots, cargo = problem
    n_slots = length(slots)
    cargo_types = cargo.cargo_types
    n_positions = length(vessel.frame_positions)
    n_cargo = length(cargo)
    n_deck = length(vessel.decks)
    n_ballast_tanks = length(vessel.ballast_tanks)

    slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)
    ρ = 1.025
    model = Model(Gurobi.Optimizer)
    # number should match number of cores used at HPC
    set_optimizer_attribute(model, "Threads", 4)
    @variable(model, weight[1:n_slots] >= 0)   # Weight at each slot
    @variable(model, cs[1:n_cargo, 1:n_slots], Bin)  # Assignment variables
    @variable(model, cargo_slack[1:n_cargo], Bin) # 1 if cargo is assigned a slot
    # Ensure placement is the same
    @constraint(model, [c = 1:n_cargo, s = 1:n_slots],
        cs[c, s] == cs1[c, s])

    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
    haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
    cost = [CSC for c in cargo]
    cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
    area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]

    # Penalty for violating constraints
    M = 100000 # Should be determined more precisely at some point
    #M = sum(cost) +1 
    # Slot weight calculation. Constraint (27)
    @constraint(model, [s = 1:n_slots],
        weight[s] == sum(cargo[c].weight * cs[c, s] for c ∈ 1:n_cargo))

    # Expression for cargo_slack
    #@expression(model,cargo_slack[c = 1:n_cargo],
    #	sum(cs[c, s] for s ∈ 1:n_slots))
    @constraint(model, [c = 1:n_cargo],
        sum(cs[c, s] for s ∈ 1:n_slots) == cargo_slack[c]
    )
    # Rasmus: Function creates a list of slot ids which are not equal to the cargo type id input "t"
    invalid_slots(t::Int) = [slot.id for slot in filter(x -> x.cargo_type_id != t, slots)]
    # Rasmus: Sets all slots which are not compatible with the cargo type to 0
    @constraint(model,
        [t in cargo_types,
            i in [cargo.id for cargo in filter(x -> x.cargo_type_id == t, cargo)]],
        sum(cs[i, j] for j in invalid_slots(t)) == 0)
    # Rasmus: Overlapping slots 
    overlapping_indicies = [Tuple(ix) for ix in findall(slots.overlap_matrix)]
    # Rasmus: Can only use one of two slots if they overlap. Constraint (26)
    @constraint(model, [(x, y) in overlapping_indicies], sum(cs[:, x]) + sum(cs[:, y]) <= 1)

    # Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
    @expression(model, pos_weight_cargo[p=1:n_positions],
        # Cargo weights using precalculated proportions
        sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots))
    # # # Deck weight limit. Constraint (29)
    @variable(model, slack_deck[1:n_deck] >= 0)
    @constraint(model, [d = 1:n_deck],
        sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit + slack_deck[d]
    )

    # Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
    @expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
    @expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
    # Rasmus: Don't understand why is not using weight[s]
    @expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

    # Slack variables:
    # Center of gravity
    @variable(model, slack_Vmax >= 0)
    @variable(model, slack_Vmin >= 0)
    @variable(model, slack_Tmin >= 0)
    @variable(model, slack_Tmax >= 0)
    @variable(model, slack_Lmin >= 0)
    @variable(model, slack_Lmax >= 0)
    # Stress and Bending
    @variable(model, slack_shear1[1:n_positions] >= 0)
    @variable(model, slack_shear2[1:n_positions] >= 0)
    stress_limits = vessel.stress_limits
    @variable(model, slack_shearMin[1:length(stress_limits)] >= 0)
    @variable(model, slack_shearMax[1:length(stress_limits)] >= 0)
    @variable(model, slack_bendingMax[1:length(stress_limits)] >= 0)
    # Ballast tanks
    @variable(model, slack_ballast_tanks[1:n_ballast_tanks] >= 0)

    add_stability_slack!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo,
        slack_Vmax, slack_Vmin, slack_Tmin, slack_Tmax, slack_Lmin, slack_Lmax,
        slack_shear1, slack_shear2, slack_shearMin, slack_shearMax, slack_bendingMax, slack_ballast_tanks)

    ballast_volume = model[:ballast_volume]

    @objective(model, Min,
        sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
        -
        sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
        +
        M * (slack_Lmax + slack_Lmin + slack_Tmax + slack_Tmin + slack_Vmax + slack_Vmin
             + sum(slack_shear1) + sum(slack_shear2) + sum(slack_shearMax) + sum(slack_shearMin) + sum(slack_bendingMax)
             + sum(slack_deck) + sum(slack_ballast_tanks))
    )
    return model
end

function model_with_slack_get_slac(sol,pro)
    model = second_stage_model_slack(sol.cs,pro)
    set_silent(model)
    set_time_limit_sec(model, 60*5) # 5 minutes
    optimize!(model)
    slack = [value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            #value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax]),
             value.(model[:slack_ballast_tanks])
            ] 
    return slack
end

# No ballast water allowed, slacked model 
function No_ballast_slacked(cs1,problem::StowageProblem)

    @unpack vessel, slots, cargo = problem
    n_slots = length(slots)
    cargo_types = cargo.cargo_types
    n_positions = length(vessel.frame_positions)
    n_cargo = length(cargo)
    n_deck = length(vessel.decks)
    n_ballast_tanks = length(vessel.ballast_tanks)

    slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)
    ρ = 1.025
    model = Model(Gurobi.Optimizer)
    # number should match number of cores used at HPC
    set_optimizer_attribute(model, "Threads", 4)
    @variable(model, weight[1:n_slots] >= 0)   # Weight at each slot
    @variable(model, cs[1:n_cargo, 1:n_slots], Bin)  # Assignment variables
    @variable(model, cargo_slack[1:n_cargo], Bin) # 1 if cargo is assigned a slot
    # Ensure placement is the same
    @constraint(model, [c = 1:n_cargo, s = 1:n_slots],
        cs[c, s] == cs1[c, s])

    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
    haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
    cost = [CSC for c in cargo]
    cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
    area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]

    # Penalty for violating constraints
    M = 100000 # Should be determined more precisely at some point
    #M = sum(cost) +1 
    # Slot weight calculation. Constraint (27)
    @constraint(model, [s = 1:n_slots],
        weight[s] == sum(cargo[c].weight * cs[c, s] for c ∈ 1:n_cargo))

    # Expression for cargo_slack
    #@expression(model,cargo_slack[c = 1:n_cargo],
    #	sum(cs[c, s] for s ∈ 1:n_slots))
    @constraint(model, [c = 1:n_cargo],
        sum(cs[c, s] for s ∈ 1:n_slots) == cargo_slack[c]
    )
    # Rasmus: Function creates a list of slot ids which are not equal to the cargo type id input "t"
    invalid_slots(t::Int) = [slot.id for slot in filter(x -> x.cargo_type_id != t, slots)]
    # Rasmus: Sets all slots which are not compatible with the cargo type to 0
    @constraint(model,
        [t in cargo_types,
            i in [cargo.id for cargo in filter(x -> x.cargo_type_id == t, cargo)]],
        sum(cs[i, j] for j in invalid_slots(t)) == 0)
    # Rasmus: Overlapping slots 
    overlapping_indicies = [Tuple(ix) for ix in findall(slots.overlap_matrix)]
    # Rasmus: Can only use one of two slots if they overlap. Constraint (26)
    @constraint(model, [(x, y) in overlapping_indicies], sum(cs[:, x]) + sum(cs[:, y]) <= 1)

    # Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
    @expression(model, pos_weight_cargo[p=1:n_positions],
        # Cargo weights using precalculated proportions
        sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots))
    # # # Deck weight limit. Constraint (29)
    @variable(model, slack_deck[1:n_deck] >= 0)
    @constraint(model, [d = 1:n_deck],
        sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit + slack_deck[d]
    )

    # Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
    @expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
    @expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
    # Rasmus: Don't understand why is not using weight[s]
    @expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

    # Slack variables:
    # Center of gravity
    @variable(model, slack_Vmax >= 0)
    @variable(model, slack_Vmin >= 0)
    @variable(model, slack_Tmin >= 0)
    @variable(model, slack_Tmax >= 0)
    @variable(model, slack_Lmin >= 0)
    @variable(model, slack_Lmax >= 0)
    # Stress and Bending
    @variable(model, slack_shear1[1:n_positions] >= 0)
    @variable(model, slack_shear2[1:n_positions] >= 0)
    stress_limits = vessel.stress_limits
    @variable(model, slack_shearMin[1:length(stress_limits)] >= 0)
    @variable(model, slack_shearMax[1:length(stress_limits)] >= 0)
    @variable(model, slack_bendingMax[1:length(stress_limits)] >= 0)
    # Ballast tanks
    @variable(model, slack_ballast_tanks[1:n_ballast_tanks] >= 0)
    # No wanting any ballast water or slack - New
    @constraint(model, [t = 1:n_ballast_tanks],
        slack_ballast_tanks[t] == 0)

    add_stability_slack_no_water!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo,
        slack_Vmax, slack_Vmin, slack_Tmin, slack_Tmax, slack_Lmin, slack_Lmax,
        slack_shear1, slack_shear2, slack_shearMin, slack_shearMax, slack_bendingMax, slack_ballast_tanks)

    ballast_volume = model[:ballast_volume]
    # New extra slack variables
    slack_dis_min = model[:slack_dis_min]
    slack_dis_max = model[:slack_dis_max]
    slack_bending_pos = model[:slack_bending_pos]
    slack_bending_neg = model[:slack_bending_neg]

    @objective(model, Min,
        sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
        -
        sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
        +
        M * (slack_Lmax + slack_Lmin + slack_Tmax + slack_Tmin + slack_Vmax + slack_Vmin
             + sum(slack_shear1) + sum(slack_shear2) + sum(slack_shearMax) + sum(slack_shearMin) + sum(slack_bendingMax)
             + sum(slack_deck) + sum(slack_ballast_tanks)+slack_dis_min+slack_dis_max+ sum(slack_bending_pos) + sum(slack_bending_neg))
    )
    return model
end