include("packages_and_files.jl")
# load data from solutions
#=
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
=#
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

# HPC_folder
test_instance = Finlandia_test[1]
HPC_folder = HPC_folders[1]#"Finlandia_"*test_instance*"_15_05_09"
#=
plot_folder = "Plots/Results/Finlandia_"*test_instance*"/"
# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end=#
# Load data from HPC
repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
println("Extra info about test: ", note)
#repetitions = 1 # didn't finish running
sc = length(scenarios)
n = length(n_unknown)
# Change if not Finlandia problem
problemname1, problemname2, problemname3 = "finlandia", test_instance, "hazardous"
Deterministic_problem = Array{Any}(undef, length(HPC_folders))
#load_data(problemname1, problemname2, problemname3)
Deterministic_Solution = Array{Any}(undef, length(HPC_folders))
#get_solution_deterministic("Finlandia_deterministic",
#    "Deterministic_Solution", HPC_folder)
# soluton arrays
#EVP_gen = Array{Any}(undef, repetitions, sc,n)
#EVP_gen_fitted = Array{Any}(undef, repetitions, sc,n)
Stochastic_gen = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
Stochastic_gen_fitted = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
#EVP_boot = Array{Any}(undef, repetitions, sc,n)
#EVP_boot_fitted = Array{Any}(undef, repetitions, sc,n)
Stochastic_boot = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
Stochastic_boot_fitted = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))

boot_fitted_slacked = Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))
gen_fitted_slacked = Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))
slack_boot= Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))
slack_gen = Array{Any}(nothing, repetitions, sc, n, length(HPC_folders))

Stochastic_problem_gen = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
Stochastic_problem_boot = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))

for l in 1:length(HPC_folders)
        test_instance = Finlandia_test[l]
        HPC_folder = HPC_folders[l]
        repetitions, scenarios, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
        println("Extra info about test: ", note)
        #repetitions = 1 # didn't finish running
        #sc = length(scenarios)
        n = length(n_unknown)
        Deterministic_Solution[l] = get_solution_deterministic("Finlandia_deterministic",
            "Deterministic_Solution", HPC_folder)
        Deterministic_problem[l] = load_data(problemname1, test_instance, problemname3)
    for i in 1:repetitions
        for j in 1:sc
            for k in 1:n
                # Solutions
                # EVP
                #foldername = "EVP_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
                #filename = "EVP_Solution"
                #EVP_gen[i,j,k] = get_solution_deterministic(foldername,
                #filename,HPC_folder)
                #EVP_gen_fitted[i,j,k] = get_solution_deterministic(foldername,
                #"Fitted_Solution",HPC_folder)
                #foldername = "EVP_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
                #filename = "EVP_Solution"
                #EVP_boot[i,j,k] = get_solution_deterministic(foldername,
                #filename,HPC_folder)
                #EVP_boot_fitted[i,j,k] = get_solution_deterministic(foldername,
                #"Fitted_Solution",HPC_folder)
                # Stochastic
                foldername = "Stochastic_random_sampling_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
                filename = "Stochastic_Problem"
                # problem
                Stochastic_problem_gen[i,j,k,l] = get_stochastic_problem(foldername,
                    filename,HPC_folders[l],problemname1,problemname2,problemname3)
                # solution
                filename = "Stochastic_Solution"
                Stochastic_gen[i, j, k,l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
                Stochastic_gen_fitted[i, j, k,l] = get_solution_deterministic(foldername,
                    "Fitted_Solution", HPC_folder)
                # Load slacked
                if Stochastic_gen_fitted[i, j, k,l].status != "OPTIMAL"
                    gen_fitted_slacked[i,j,k,l] = get_solution_deterministic(foldername,
                    "Fitted_Solution_slacked", HPC_folder)
                    slack_gen[i,j,k,l] = get_slack(foldername, "Fitted_Solution",HPC_folder)
        
                end
                # Stochastic Bootstrap
                foldername = "Stochastic_Bootstrap1_rep$(i)_sc$(scenarios[j])_unknown$(n_unknown[k])_time$(time_limit)"
                filename = "Stochastic_Problem"
                # problem
                Stochastic_problem_boot[i,j,k,l] = get_stochastic_problem(foldername,
                    filename,HPC_folders[l],problemname1,problemname2,problemname3)
                # solution
                filename = "Stochastic_Solution"
                Stochastic_boot[i, j, k,l] = get_solution_stochastic(foldername,
                    filename, HPC_folder)
                Stochastic_boot_fitted[i, j, k,l] = get_solution_deterministic(foldername,
                    "Fitted_Solution", HPC_folder)
                # Load slacked
                if Stochastic_boot_fitted[i, j, k,l].status != "OPTIMAL"
                    boot_fitted_slacked[i,j,k,l] = get_solution_deterministic(foldername,
                    "Fitted_Solution_slacked", HPC_folder)
                    slack_boot[i,j,k,l] = get_slack(foldername, "Fitted_Solution",HPC_folder)
                end
            end
        end
    end
end
# plot 
plot_folder = "Plots/Results/Finlandia_plots/"
gaps_boot = zeros(repetitions, sc,length(HPC_folders))
gaps_gen = zeros(repetitions, sc,length(HPC_folders))
cargo_loaded_boot = zeros(repetitions, sc,length(HPC_folders))
cargo_loaded_gen = zeros(repetitions, sc,length(HPC_folders))
ballast_water_boot = zeros(repetitions, sc,length(HPC_folders))
ballast_water_gen = zeros(repetitions, sc,length(HPC_folders))
for l in 1:length(HPC_folders)
    for i in 1:repetitions
        for j in 1:sc
            if Stochastic_boot_fitted[i, j, 1,l].status != "OPTIMAL"
                println("Boot - status: ", Stochastic_boot_fitted[i, j, 1,l].status)
                println("test: $(l), rep: $(i), sc: $(j)")
            end
            if Stochastic_gen_fitted[i,j,1,l].status != "OPTIMAL"
                println("gen - status: ", Stochastic_gen_fitted[i, j, 1,l].status)
                println("test: $(l), rep: $(i), sc: $(j)")
            end

            gaps_boot[i, j,l] = round(Stochastic_boot[i, j, 1,l].gap*100,digits = 3)
            gaps_gen[i, j,l] = round(Stochastic_gen[i, j, 1,l].gap*100,digits = 3)
            cargo_loaded_boot[i, j,l] = Int(Stochastic_boot[i, j, 1,l].n_cargo_loaded)
            cargo_loaded_gen[i, j,l] = Int(Stochastic_gen[i, j, 1,l].n_cargo_loaded)
            ballast_water_boot[i, j,l] = Stochastic_boot_fitted[i, j, 1,l].ballast_weight
            ballast_water_gen[i, j,l] = Stochastic_gen_fitted[i, j, 1,l].ballast_weight
        end
    end
end
means_gen = mean(gaps_gen,dims=1)
light_rows = [1,2,5,6]
heavy_rows = [3,4,7,8]
means_gen_light = dropdims(means_gen[:,:,light_rows];dims=1)
mean_heavy_gen = dropdims(means_gen[:,:,heavy_rows];dims=1)
mean_diff_gen = means_gen_light .- mean_heavy_gen

pvalues_gen = zeros(sc,length(light_rows))
pvalues_boot = zeros(sc,length(light_rows))
for i in 1:length(light_rows)
    for j in 1:sc
        x1 = gaps_gen[:,j,light_rows[i]]
        x2 = gaps_gen[:,j,heavy_rows[i]]
        t_test = EqualVarianceTTest(x1, x2)
        pvalues_gen[j,i] = round(pvalue(t_test),digits=4)
        if pvalues_gen[j,i] < 0.05
            println("Instance $(light_rows[i]) and $(heavy_rows[i]) are significantly different for scenario $(j) with p-value: ", pvalues_gen[j,i])
        end
        x1 = gaps_boot[:,j,light_rows[i]]
        x2 = gaps_boot[:,j,heavy_rows[i]]
        t_test = EqualVarianceTTest(x1, x2)
        pvalues_boot[j,i] = round(pvalue(t_test),digits=4)
        if pvalues_boot[j,i] < 0.05
            println("Instance $(light_rows[i]) and $(heavy_rows[i]) are significantly different for scenario $(j) with p-value: ", pvalues_boot[j,i])
        end
    end
end 
pretty_table(pvalues_gen)
pretty_table(pvalues_boot)

# One test is fucked - running same problem through HPC again to see if it happens again
getfield.(Stochastic_boot[:,5,1,8], :n_cargo_loaded)

for i in 1:length(HPC_folders)
    #println("##############")
    println("Instance: $(i)")
    println("Det gap: ", Deterministic_Solution[i].gap*100)
    println("Det Status: ", Deterministic_Solution[i].status)
    #println(round(Deterministic_Solution[i].ballast_weight,digits=2))
    #println("Gap: ",Deterministic_Solution[i].gap)
    #println("Cargo loaded: ", Deterministic_Solution[i].n_cargo_loaded,"/",
    #Deterministic_Solution[i].n_cargo_total)
    for j in 1:length(scenarios)
        #println("Gap, sc: $(j): ", mean(gaps_gen[:,j,i]))
        #println("SD, sc: $(j): ", round(std(gaps_boot[:,j,i]),digits=2))
        #println("Cargo loaded, sc: $(j): ", mean(cargo_loaded_boot[:,j,i]))
        #println("Ballast water, sc: $(j): ", round(mean(ballast_water_boot[:,j,i]),digits=2))
        #println("Feasible, sc: $(j): ", count(x-> x== nothing, boot_fitted_slacked[:,j,1,i]),"/",repetitions)
    end
end
count(x-> x!= nothing, boot_fitted_slacked[:,5,1,8])
# Finding ballast water for tests where some instances were infeasible
println(ballast_water_boot[:,1,7])
# instance 7 boot
println(getfield.(Stochastic_boot_fitted[:,1,1,7],:status))
round(mean(filter(x->x>0, ballast_water_boot[:,1,7])),digits=2)
Stochastic_boot_fitted[5, 1, 1,7].status
# instance 8 boot - 40 scenarios
println(getfield.(Stochastic_boot_fitted[:,4,1,8],:status))
println(ballast_water_boot[:,4,8])
round(mean(filter(x->x>0, ballast_water_boot[:,4,8])),digits=2)
# instance 8 boot - 50 scenarios
println(getfield.(Stochastic_boot_fitted[:,5,1,8],:status))
println(ballast_water_boot[:,5,8])
round(mean(filter(x->x>0, ballast_water_boot[:,5,8])),digits=2)


Stochastic_boot[9,5,1,8].cargo
Stochastic_problem_boot[9,5,1,8].cargo.items.total_weight
println(getfield.(Stochastic_gen[:,5,1,8],:n_cargo_loaded))
println(getfield.(Stochastic_boot[:,5,1,8],:n_cargo_loaded))





# Sizes of the models
n_var = []
n_const = []
m_size = []
for i in 1:sc
    push!(n_var, Stochastic_gen[1,i,1,2].n_variables)
    push!(n_const, Stochastic_gen[1,i,1,2].n_constraints)
    push!(m_size, Stochastic_gen[1,i,1,2].model_size)
end
p1= plot(scenarios,n_var, title = "Number of variables in instance 2",
xlabel = "Scenarios", ylabel = "Number of variables", label = "")
savefig(p1, "Plots/Data/"*"NumberOfVariables_Instance2.png")
p2 = plot(scenarios,n_const, title = "Number of constraints in instance 2",
xlabel = "Scenarios", ylabel = "Number of constraints", label = "")
savefig(p2, "Plots/Data/"*"NumberOfConstraints_Instance2.png")
p3= plot(scenarios,m_size, title = "Model size in instance 2",
xlabel = "Scenarios", ylabel = "Model size (Bytes)", label = "")#,yscale=:log10)
savefig(p3, "Plots/Data/"*"ModelSize_Instance2.png")

# Infeasible models
slack_string = ["deck", "Vmax", "Vmin", "Tmin", "Tmax", "Lmin", "Lmax",
"shear1", "shear2",
"shearMin", "shearMax", "bendingMax", "ballast_tanks"]
slack_count_gen = zeros(length(slack_string))
slack_count_boot = zeros(length(slack_string))
deck_gen = zeros(3)
deck_boot = zeros(3)
deck_1_violations_gen = zeros(repetitions,sc,length(HPC_folders))
deck_1_violations_boot = zeros(repetitions,sc,length(HPC_folders))
for i in 1:repetitions
    for j in 1:sc
        for l in 1:length(HPC_folders)
            if !isnothing(boot_fitted_slacked[i,j,1,l])
                println("#####################")
                println("Repetition: $(i), Scenarios: $(j), Instance: $(l), Boot fitted slacked is not nothing")
                for k in 1:length(slack_boot[i,j,1,l])
                    if sum(slack_boot[i,j,1,l][k]) > 0 
                        println("slack variable different than 0: ", slack_string[k])
                        println("Slack: ", slack_boot[i,j,1,l][k])
                        slack_count_boot[k] += 1
                        if k==1
                            deck_boot[k] += 1
                            deck_1_violations_boot[i,j,l] = sum(slack_boot[i,j,1,l][k])
                        end
                    end
                end
            end
            if !isnothing(gen_fitted_slacked[i,j,1,l])
                println("#####################")
                println("Repetition: $(i), Scenarios: $(j), Instance: $(l), gen fitted slacked is not nothing")
                for k in 1:length(slack_gen[i,j,1,l])
                    if sum(slack_gen[i,j,1,l][k]) > 0 
                        println("slack variable different than 0: ", slack_string[k])
                        println("Slack: ", slack_gen[i,j,1,l][k])
                        slack_count_gen[k] += 1
                        if k==1
                            deck_gen[k] += 1
                            deck_1_violations_gen[i,j,l] = sum(slack_gen[i,j,1,l][k])
                        end
                    end
                end
            end
        end
    end
end
# slack variable order:
# deck, Vmax, Vmin, Tmin, Tmax, Lmin, Lmax, shear1, shear2, 
#shearMin, shearMax, bendingMax, ballast_tanks
println(slack_count_gen)
println(slack_count_boot)
println(deck_gen)
println(deck_boot)

mean(deck_1_violations_gen[:,1,8])
mean(deck_1_violations_boot[:,1,8])
mean(deck_1_violations_gen[:,2,8])
mean(deck_1_violations_boot[:,2,8])
mean(deck_1_violations_gen[:,3,8])
mean(deck_1_violations_boot[:,3,8])
#for i in 1:repetitions

println(getfield.(Stochastic_boot[:,5,1,8],:n_cargo_loaded))
println(getfield.(Stochastic_boot_fitted[:,5,1,8],:status))
Stochastic_boot_fitted[3,5,1,8].ballast_weight
Stochastic_boot_fitted[9,5,1,8].ballast_weight
Stochastic_boot_fitted[10,5,1,8].ballast_weight


slots_for_trucks = filter(x->x.cargo_type_id == 1, Deterministic_problem[8].slots)
slots_for_trucks_deck1 = filter(x->x.deck_id == 1, slots_for_trucks)
length(slots_for_trucks)
length(slots_for_trucks_deck1)
trucks_placed = filter(x->x.cargo_type_id==1, Deterministic_Solution[8].cargo)
trucks_placed_on_dekc1 = filter(x->x.deck==1, trucks_placed)
println("Trucks placed: ", length(trucks_placed),"/",
length(filter(x->x.cargo_type_id == 1, Deterministic_problem[8].cargo)))
println("Trucks placed on deck 1: ", length(trucks_placed_on_dekc1),"/",
length(slots_for_trucks_deck1))
mean_weight_truck = mean(c.weight for c in filter(x->x.cargo_type_id == 1, Deterministic_problem[8].cargo))
println("Mean weight of truck: ", mean_weight_truck)
println("Mean weight trucks times slots for trucks on deck 1 / weight limit deck 1: ",
mean_weight_truck*length(slots_for_trucks_deck1),"/", Deterministic_problem[8].vessel.decks[1].weight_limit)
mean_weight_truck_weight_interval = (40-5)/2
println("Mean weight trucks times slots for trucks on deck 1 / weight limit deck 1: ",
mean_weight_truck_weight_interval*length(slots_for_trucks_deck1),"/", Deterministic_problem[8].vessel.decks[1].weight_limit)


#=
bplots_gen = []
bplots_boot = []
for i in 1:length(HPC_folders)
    p = boxplot(gaps_gen[:,1,i], label = "No. Scenarios 10", 
    ylabel="Gap", title="Gap for instance $(i),\n SGM. 1")
    boxplot!(p,gaps_gen[:,2,i], label = "No. Scenarios 20")
    boxplot!(p,gaps_gen[:,3,i], label = "No. Scenarios 30")
    boxplot!(p,gaps_gen[:,4,i], label = "No. Scenarios 40")
    boxplot!(p,gaps_gen[:,5,i], label = "No. Scenarios 50")
    push!(bplots_gen, p)
    p = boxplot(gaps_boot[:,1,i], label = "No. Scenarios 10", 
    ylabel="Gap", title="Gap for instance $(i),\n SGM. 2")
    boxplot!(p,gaps_boot[:,2,i], label = "No. Scenarios 20")
    boxplot!(p,gaps_boot[:,3,i], label = "No. Scenarios 30")
    boxplot!(p,gaps_boot[:,4,i], label = "No. Scenarios 40")
    boxplot!(p,gaps_boot[:,5,i], label = "No. Scenarios 50")
    push!(bplots_boot, p)   
end
display(bplots_gen[2])
display(bplots_boot[2])
=#

# Scenarios
# Comparing instance 1 and 2 in SGM 1 and 2
plot_folder_scenarios = "Plots/Data/Scenarios/"
foldername_boot = "Stochastic_Bootstrap1_rep$(1)_sc$(scenarios[5])_unknown$(70)_time$(time_limit)"
foldername_gen = "Stochastic_random_sampling_rep$(1)_sc$(scenarios[5])_unknown$(70)_time$(time_limit)"
filename = "Stochastic_Problem"
problemname1, problemname2, problemname3 = "finlandia", Finlandia_test[5], "hazardous"
gen_problem = get_stochastic_problem(foldername_gen,filename,
            HPC_folders[5],problemname1,problemname2,problemname3)
boot_problem = get_stochastic_problem(foldername_boot,filename,
            HPC_folders[5],problemname1,problemname2,problemname3)
weight_diff_gen = gen_problem.cargo.items.total_weight .- ones(gen_problem.scenarios)*Deterministic_problem[5].cargo.total_weight
weight_diff_boot = boot_problem.cargo.items.total_weight .- ones(boot_problem.scenarios)*Deterministic_problem[5].cargo.total_weight
p1 = plot(weight_diff_gen, title = "Total cargo weight difference in each scenario",
xlabel = "Scenario", ylabel = "Weight difference (t)", label = "")
savefig(p1, plot_folder_scenarios*"WeightDiff_Gen_instance_5.png")
p2 = plot(weight_diff_boot, title = "Total cargo weight difference in each scenario",
xlabel = "Scenario", ylabel = "Weight difference (t)", label = "")
savefig(p2, plot_folder_scenarios*"WeightDiff_Boot_instance_5.png")
std_5 = [std(weight_diff_gen),std(weight_diff_boot)]

cargo_1_weight_gen = [gen_problem.cargo.items[i].items[25].weight for i in 1:gen_problem.scenarios]
cargo_1_weight_boot = [boot_problem.cargo.items[i].items[25].weight for i in 1:boot_problem.scenarios]
p3 = scatter(cargo_1_weight_gen, title = "Cargo 25's weight in each scenario",
xlabel = "Scenario", ylabel = "Weight (t)", label = "Scenario weight")
plot!(p3, ones(length(cargo_1_weight_gen))*Deterministic_problem[5].cargo.items[25].weight,
label = "Actual weight")
savefig(p3, plot_folder_scenarios*"Cargo1Weight_Gen_instance_5.png")
p4 = scatter(cargo_1_weight_boot, title = "Cargo 25's weight in each scenario",
xlabel = "Scenario", ylabel = "Weight (t)", label = "")
plot!(p4, ones(length(cargo_1_weight_boot))*Deterministic_problem[5].cargo.items[25].weight,
label = "Actual weight")
savefig(p4, plot_folder_scenarios*"Cargo1Weight_Boot_instance_5.png")

# Heavier instance
foldername_boot = "Stochastic_Bootstrap1_rep$(1)_sc$(scenarios[5])_unknown$(70)_time$(time_limit)"
foldername_gen = "Stochastic_random_sampling_rep$(1)_sc$(scenarios[5])_unknown$(70)_time$(time_limit)"
filename = "Stochastic_Problem"
problemname1, problemname2, problemname3 = "finlandia", Finlandia_test[7], "hazardous"
gen_problem = get_stochastic_problem(foldername_gen,filename,
            HPC_folders[7],problemname1,problemname2,problemname3)
boot_problem = get_stochastic_problem(foldername_boot,filename,
            HPC_folders[7],problemname1,problemname2,problemname3)
weight_diff_gen = gen_problem.cargo.items.total_weight .- ones(gen_problem.scenarios)*Deterministic_problem[7].cargo.total_weight
weight_diff_boot = boot_problem.cargo.items.total_weight .- ones(boot_problem.scenarios)*Deterministic_problem[7].cargo.total_weight
p1 = plot(weight_diff_gen, title = "Total cargo weight difference in each scenario",
xlabel = "Scenario", ylabel = "Weight difference (t)", label = "")
savefig(p1, plot_folder_scenarios*"WeightDiff_Gen_instance_7.png")
p2 = plot(weight_diff_boot, title = "Total cargo weight difference in each scenario",
xlabel = "Scenario", ylabel = "Weight difference (t)", label = "")
savefig(p2, plot_folder_scenarios*"WeightDiff_Boot_instance_7.png")
std_7 = [std(weight_diff_gen),std(weight_diff_boot)]


cargo_1_weight_gen = [gen_problem.cargo.items[i].items[25].weight for i in 1:gen_problem.scenarios]
cargo_1_weight_boot = [boot_problem.cargo.items[i].items[25].weight for i in 1:boot_problem.scenarios]
p3 = scatter(cargo_1_weight_gen, title = "Cargo 25's weight in each scenario",
xlabel = "Scenario", ylabel = "Weight (t)", label = "Scenario weight")
plot!(p3, ones(length(cargo_1_weight_gen))*Deterministic_problem[7].cargo.items[25].weight,
label = "Actual weight")
savefig(p3, plot_folder_scenarios*"Cargo1Weight_Gen_instance_7.png")
p4 = scatter(cargo_1_weight_boot, title = "Cargo 25's weight in each scenario",
xlabel = "Scenario", ylabel = "Weight (t)", label = "")
plot!(p4, ones(length(cargo_1_weight_boot))*Deterministic_problem[7].cargo.items[25].weight,
label = "Actual weight")
savefig(p4, plot_folder_scenarios*"Cargo1Weight_Boot_instance_7.png")

std_5
std_7

for i in 1:gen_problem.scenarios
    temp = [c.weight for c in( filter(x->x.cargo_type_id == 1, gen_problem.cargo.items[i]))]
    println("mean: ",mean(temp))
    println("max: ",maximum(temp))
    println("min: ",minimum(temp))
    println(argmax(temp))
end
temp = [c.weight for c in( filter(x->x.cargo_type_id == 1, Deterministic_problem[7].cargo))]
a = argmax(temp)
minimum(temp)
gen_problem.cargo.items[1].items[32].cargo_type_id
[gen_problem.cargo.items[i].items[32].cargo for i in 1:gen_problem.scenarios]

for i in 1:length(HPC_folders)
    # print cargo Weight
    #println("Cargo weight in problem $(i): ",round(Deterministic_problem[i].cargo.total_weight,digits=2))
    #println("Avg. cargo weight - gen, instance $(i): ",
    #round(mean(Stochastic_problem_gen[1,end,1,i].cargo.items.total_weight),digits=2))
    println("Avg. cargo weight - boot, instance $(i): ",
    round(mean(Stochastic_problem_boot[1,end,1,i].cargo.items.total_weight),digits=2))
end

#Stochastic_problem_gen = Array{Any}(undef, repetitions, sc, n, length(HPC_folders))
















#############################
# Old tests
for l in 1:length(HPC_folders)
    println("#########################")
    println("Problem:", Finlandia_test[l])
    println("Problem idx: $(l)")
    #println("Gaps - boot:")
    pretty_table(gaps_boot[:,:,l], 
    header = ["sc: $(scenarios[i])" for i in 1:sc], 
    row_labels = ["rep: $(i)" for i in 1:repetitions],title = "Boot - Gaps")
    #println("Gaps - gen:")
    pretty_table(gaps_gen[:,:,l], 
    header = ["sc: $(scenarios[i])" for i in 1:sc], 
    row_labels = ["rep: $(i)" for i in 1:repetitions],title = "Gen - Gaps")
    #println("Cargo loaded - boot:")
    pretty_table(cargo_loaded_boot[:,:,l], 
    header = ["sc: $(scenarios[i])" for i in 1:sc], 
    row_labels = ["rep: $(i)" for i in 1:repetitions], title="Boot - Cargo loaded. Det_sol: $(Deterministic_Solution[l].n_cargo_loaded)/$(Deterministic_Solution[l].n_cargo_total)")
    #println("Cargo loaded - gen:")
    pretty_table(cargo_loaded_gen[:,:,l], 
    header = ["sc: $(scenarios[i])" for i in 1:sc], 
    row_labels = ["rep: $(i)" for i in 1:repetitions], title="Gen - Cargo loaded. Det_sol: $(Deterministic_Solution[l].n_cargo_loaded)/$(Deterministic_Solution[l].n_cargo_total)")
end


#=
include("src/utils/SaveData.jl")
test = get_stochastic_problem("Stochastic_Bootstrap1_rep$(1)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)",
"Stochastic_Problem",HPC_folder,problemname1,problemname2,problemname3)
cargo = test.cargo
test3 = create_model_stochastic_cargo_fraction(test,0.9)
set_silent(test2) # removes terminal output
set_time_limit_sec(test3, 5*60) # 5 minutes to start with
optimize!(test3)
=#

#test_Det = get_solution_deterministic("Finlandia_deterministic","Deterministic_Solution",HPC_folder)

println("##########################")

println("Cargo in problem: ", length(Deterministic_problem.cargo))
println("Deterministic solution packs: ",Deterministic_Solution.n_cargo_loaded)

println("##########################")
println("Total number of models :", repetitions*sc*n)
println("Deterministic solution, Cargo loaded: ", Deterministic_Solution.n_cargo_loaded)
for i in 1:repetitions
    for j in 1:sc
        for k in 1:n
            println("############################")
            println("Number of scnearios: ", scenarios[j])
            #println("EVP model gen, Cargo loaded: ", EVP_gen[i,j,k].n_cargo_loaded)
            #println("EVP model gen after realization, Cargo loaded: ", EVP_gen_fitted[i,j,k].n_cargo_loaded)
            println("Stochastic model gen, Cargo loaded: ", Stochastic_gen[i,j,k].n_cargo_loaded)
            println("Stochastic model gen after realization, Cargo loaded: ", Stochastic_gen_fitted[i,j,k].n_cargo_loaded)
            #println("EVP model boot, Cargo loaded: ", EVP_boot[i,j,k].n_cargo_loaded)
            #println("EVP model boot after realization, Cargo loaded: ", EVP_boot_fitted[i,j,k].n_cargo_loaded)
            println("Stochastic model boot, Cargo loaded: ", Stochastic_boot[i,j,k].n_cargo_loaded)
            println("Stochastic model boot after realization, Cargo loaded: ", Stochastic_boot_fitted[i,j,k].n_cargo_loaded)
        end
    end
end
println("##########################")
for i in 1:sc
    println("Average number of cargo loaded, scenarios = $(scenarios[i]) - gen: ",
    sum([Stochastic_gen[j,i,end].n_cargo_loaded for j in 1:repetitions])/repetitions)
    println("Average number of cargo loaded, scenarios = $(scenarios[i]) - boot: ",
    sum([Stochastic_boot[j,i,end].n_cargo_loaded for j in 1:repetitions])/repetitions)
end
println("##########################")

# Plot cargo loaded
n_models = 4
r = 1 # choose repetition number
#cargo_loaded_EVP_gen = zeros(Int64,sc,n)
#cargo_loaded_EVP_gen_fitted = zeros(Int64,sc,n)
cargo_loaded_Stochastic_gen = zeros(Int64,sc,n)
cargo_loaded_Stochastic_gen_fitted = zeros(Int64,sc,n)
#cargo_loaded_EVP_boot = zeros(Int64,sc,n)
#cargo_loaded_EVP_boot_fitted = zeros(Int64,sc,n)
cargo_loaded_Stochastic_boot = zeros(Int64,sc,n)
cargo_loaded_Stochastic_boot_fitted = zeros(Int64,sc,n)
#gap_EVP_gen = zeros(Float64,sc,n)
gap_Stochastic_gen = zeros(Float64,sc,n)
#gap_EVP_boot = zeros(Float64,sc,n)
gap_Stochastic_boot = zeros(Float64,sc,n)
for i in 1:sc
    for j in 1:n
        #cargo_loaded_EVP_gen[i,j] = Int(EVP_gen[r,i,j].n_cargo_loaded)
        #cargo_loaded_EVP_gen_fitted[i,j] = Int(EVP_gen_fitted[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_gen[i,j] = Int(Stochastic_gen[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_gen_fitted[i,j] = Int(Stochastic_gen_fitted[r,i,j].n_cargo_loaded)
        #cargo_loaded_EVP_boot[i,j] = Int(EVP_boot[r,i,j].n_cargo_loaded)
        #cargo_loaded_EVP_boot_fitted[i,j] = Int(EVP_boot_fitted[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_boot[i,j] = Int(Stochastic_boot[r,i,j].n_cargo_loaded)
        cargo_loaded_Stochastic_boot_fitted[i,j] = Int(Stochastic_boot_fitted[r,i,j].n_cargo_loaded)
        #gap_EVP_gen[i,j] = EVP_gen[r,i,j].gap
        gap_Stochastic_gen[i,j] = Stochastic_gen[r,i,j].gap
        #gap_EVP_boot[i,j] = EVP_boot[r,i,j].gap
        gap_Stochastic_boot[i,j] = Stochastic_boot[r,i,j].gap
    end
end
# Display loaded cargo
# EVP gen
linew = 1
p = plot(scenarios,cargo_loaded_Stochastic_gen[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Cargo loaded for different models,\n before realization of cargo weight", label = "Sto-Gen", linewidth =linew)
#plot!(scenarios,cargo_loaded_EVP_boot[:,end],label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,cargo_loaded_Stochastic_gen[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_boot[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.n_cargo_loaded,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"CargoLoaded_BeforeRealization.png")
p = plot(scenarios,cargo_loaded_Stochastic_gen_fitted[:,end],xlabel="Scenarios",ylabel="Cargo loaded",
title="Cargo loaded for different models,\n after realization of cargo weight", label = "Sto-Gen", linewidth =linew)
#plot!(scenarios,cargo_loaded_EVP_boot_fitted[:,end],label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,cargo_loaded_Stochastic_gen_fitted[:,end],label = "Sto-Gen", linewidth =linew)
plot!(scenarios,cargo_loaded_Stochastic_boot_fitted[:,end],label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.n_cargo_loaded,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"CargoLoaded_AfterRealization.png")

# Display Gaps
p = plot(scenarios,gap_Stochastic_gen[:,end] .*100,xlabel="Scenarios",ylabel="Gap (%)",
title="Gap for different models,\n before realization of cargo weight", label = "Sto-Gen", linewidth =linew)
#plot!(scenarios,gap_EVP_boot[:,end] .*100,label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,gap_Stochastic_gen[:,end] .* 100,label = "Sto-Gen", linewidth =linew)
plot!(scenarios,gap_Stochastic_boot[:,end] .* 100,label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.gap .*100,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"Gap_BeforeRealization.png")
# Zoomed in
p = plot(scenarios,gap_Stochastic_boot[:,end].*100,xlabel="Scenarios",ylabel="Gap (%)",
title="Gap for different models,\n before realization of cargo weight", label = "Sto-Boot", linewidth =linew)
#plot!(scenarios,gap_EVP_boot[:,end].*100,label = "EVP-Boot", linewidth =linew)
#plot!(scenarios,gap_Stochastic_gen[:,end],label = "Sto-Gen", linewidth =linew)
#plot!(scenarios,gap_Stochastic_boot[:,end].*100,label = "Sto-Boot", linewidth =linew)
plot!(scenarios,ones(length(scenarios))*Deterministic_Solution.gap .*100,label = "Deterministic", linewidth =linew, linestyle=:dash)
display(p)
savefig(plot_folder*"Gap_BeforeRealization_zoomed.png")

##########################
# Display ballast water from different models
p = plot(scenarios, ones(sc)*Deterministic_Solution.ballast_weight, marker = :utriangle, xlabel = "Scenarios", 
ylabel = "Ballast weight (t)", label = "Deterministic", title = "Ballast weight for different models",
formatter=:plain)
annotate!(scenarios[1], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))
annotate!(scenarios[2], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))
annotate!(scenarios[3], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))
annotate!(scenarios[4], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))
annotate!(scenarios[5], Deterministic_Solution.ballast_weight + 15, text(Deterministic_Solution.n_cargo_loaded, 6))

xtemp = [[],[],[],[]]
ytemp = [[],[],[],[]]
cargo_n = [[],[],[],[]]
for i in 1:sc
    if Stochastic_gen_fitted[1,i,end].status != "INFEASIBLE" # have a solution 
        push!(xtemp[1],scenarios[i])
        push!(ytemp[1],Stochastic_gen_fitted[1,i,end].ballast_weight)
        push!(cargo_n[1],Stochastic_gen_fitted[1,i,end].n_cargo_loaded)
    end
    if Stochastic_boot_fitted[1,i,end].status != "INFEASIBLE" # have a solution 
        push!(xtemp[2],scenarios[i])
        push!(ytemp[2],Stochastic_boot_fitted[1,i,end].ballast_weight)
        push!(cargo_n[2],Stochastic_boot_fitted[1,i,end].n_cargo_loaded)
    end
    #if EVP_gen_fitted[1,i,end].status != "INFEASIBLE" # have a solution 
    #    push!(xtemp[3],scenarios[i])
    #    push!(ytemp[3],EVP_gen_fitted[1,i,end].ballast_weight)
    #    push!(cargo_n[3],EVP_gen_fitted[1,i,end].n_cargo_loaded)
    #end
    #if EVP_boot_fitted[1,i,end].status != "INFEASIBLE" # have a solution 
    #    push!(xtemp[4],scenarios[i])
    #    push!(ytemp[4],EVP_boot_fitted[1,i,end].ballast_weight)
    #    push!(cargo_n[4],EVP_boot_fitted[1,i,end].n_cargo_loaded)
    #end
end
plot!(xtemp[1],ytemp[1], label = "Stochastic Gen", marker = :circle, markersize = 4)
plot!(xtemp[2],ytemp[2], label = "Stochastic Boot", marker = :circle, markersize = 4)
plot!(formatte=:plain)
#plot!(xtemp[3],ytemp[3], label = "EVP Gen", marker = :circle, markersize = 4)
#plot!(xtemp[4],ytemp[4], label = "EVP Boot", marker = :circle, markersize = 4)
for j in 1:length(xtemp)
    for i = 1:length(xtemp[j])
        annotate!(xtemp[j][i], ytemp[j][i] + 20, text(cargo_n[j][i],6))
    end
end
display(p)
savefig(plot_folder*"BallastWeight.png")

# Number of times problem was unfeasible
#EVP_gen_inf = []
#EVP_gen_fitted_inf = []
Stochastic_gen_inf = []
Stochastic_gen_fitted_inf = []
#EVP_boot_inf = []
#EVP_boot_fitted_inf = []
Stochastic_boot_inf = []
Stochastic_boot_fitted_inf = []
println("Checing problems for infeasibility")
for i in 1:sc
    for j in 1:n
        # print if not optimal
        #if EVP_gen[r,i,j].status != "OPTIMAL"
        #    push!(EVP_gen_inf,EVP_gen[r,i,j])
        #    println("##########################")
        #    println("EVP gen model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
        #    println("Status: ", EVP_gen[r,i,j].status)
        #end
        #if EVP_gen_fitted[r,i,j].status != "OPTIMAL"
        #    push!(EVP_gen_fitted_inf,EVP_gen_fitted[r,i,j])
        #    println("##########################")
        #    println("EVP gen fitted model jnfeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
        #    println("Status: ", EVP_gen_fitted[r,i,j].status)
        #end
        if Stochastic_gen[r,i,j].status != "OPTIMAL"
            push!(Stochastic_gen_inf,Stochastic_gen[r,i,j])
            println("##########################")
            println("Stochastic gen model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
            println("Status: ", Stochastic_gen[r,i,j].status)
        end
        if Stochastic_gen_fitted[r,i,j].status != "OPTIMAL"
            push!(Stochastic_gen_fitted_inf,Stochastic_gen_fitted[r,i,j])
            println("##########################")
            println("Stochastic gen fitted model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
            println("Status: ", Stochastic_gen_fitted[r,i,j].status)
        end
        #if EVP_boot[r,i,j].status != "OPTIMAL"
        #    push!(EVP_boot_inf,EVP_boot[r,i,j])
        #    println("##########################")
        #    println("EVP boot model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
        #    println("Status: ", EVP_boot[r,i,j].status)
        #end
        #if EVP_boot_fitted[r,i,j].status != "OPTIMAL"
        #    push!(EVP_boot_fitted_inf,EVP_boot_fitted[r,i,j])
        #    println("##########################")
        #    println("EVP boot fitted model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
        #    println("Status: ", EVP_boot_fitted[r,i,j].status)
        #end
        if Stochastic_boot[r,i,j].status != "OPTIMAL"
            push!(Stochastic_boot_inf,Stochastic_boot[r,i,j])
            println("##########################")
            println("Stochastic boot model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
            println("Status: ", Stochastic_boot[r,i,j].status)
        end
        if Stochastic_boot_fitted[r,i,j].status != "OPTIMAL"
            push!(Stochastic_boot_fitted_inf,Stochastic_boot_fitted[r,i,j])
            println("##########################")
            println("Stochastic boot fitted model infeasible for repetition $(r), scenario $(scenarios[(i)]), unknown cargo $(j)")
            println("Status: ", Stochastic_boot_fitted[r,i,j].status)
        end
    end
end
println("##########################")
println("Number of different parameters: ", sc*n)
#println("Stochastic models + EVP models: ", length(EVP_gen)+length(EVP_boot)+length(Stochastic_gen)+length(Stochastic_boot))
println("##########################")
# EVP models - normally none of them should be infeasible
#=
for i in 1:length(EVP_gen_inf)
    println("EVP-gen Status: ", EVP_gen_inf[i].status)
    println("Index: ", i)
end
println("##########################")
for i in 1:length(EVP_boot_inf)
    println("EVP-boot Status: ", EVP_boot_inf[i].status)
    println("Index: ", i)
end
=#
# Stochastic models 
println("##########################")
println("Stochastic gen model:")
for i in 1:length(Stochastic_gen_inf)
    println("Stochastic gen model infeasible for number of scenarios: $(length(Stochastic_gen_inf[i].forces))")
    println("Status: ", Stochastic_gen_inf[i].status)
    if Stochastic_gen_inf[i].status == "TIME_LIMIT"
        println("Cargo loaded: ", Stochastic_gen_inf[i].n_cargo_loaded)
    else # if infeasible
        for j in 1:lengt(h(Stochastic_gen_inf[i].forces))
            println("Weight in scenarios", Stochastic_gen_inf[i].forces[j].cargo_weight)
        end
    end
end
println("##########################")
println("Stochastic Boot model:")
for i in 1:length(Stochastic_boot_inf)
    println("Stochastic boot model infeasible for number of scenarios: $(length(Stochastic_boot_inf[i].forces))")
    println("Status: ", Stochastic_boot_inf[i].status)
    if Stochastic_boot_inf[i].status == "TIME_LIMIT"
        println("Cargo loaded: ", Stochastic_boot_inf[i].n_cargo_loaded)
    else # if infeasible
        for j in 1:lengt(h(Stochastic_boot_inf[i].forces))
            println("Weight in scenarios", Stochastic_boot_inf[i].forces[j].cargo_weight)
        end
    end
end
# After the realization: The recourse Model
println("##########################")
# EVP models - normally none of them should be infeasible
#=
println("EVP-fitted gen model:")
for i in 1:length(EVP_gen_fitted_inf)
    println("Status: ", EVP_gen_fitted_inf[i].status)
end
println("##########################")
println("EVP-fitted boot model:")
for i in 1:length(EVP_boot_fitted_inf)
    println("Status: ", EVP_boot_fitted_inf[i].status)
    if EVP_boot_fitted_inf[i].status != "TIME_LIMIT"
        # TODO Do something to check why this is the case
    end
end
=#
# Stochastic models 
println("##########################")
println("Stochastic-fitted gen model:")
for i in 1:length(Stochastic_gen_fitted_inf)
    println("Stochastic gen model infeasible.")
    println("Status: ", Stochastic_gen_fitted_inf[i].status)
    if Stochastic_gen_fitted_inf[i].status != "TIME_LIMIT" # infeasible
        # TODO Do something to check why this is the case
    end
end
println("##########################")
println("Stochastic-fitted boot model:")
for i in 1:length(Stochastic_boot_fitted_inf)
    println("Stochastic boot model infeasible.")
    println("Status: ", Stochastic_boot_fitted_inf[i].status)
    if Stochastic_boot_fitted_inf[i].status != "TIME_LIMIT"
        # TODO Do something to check why this is the case
    end
end

##################################3
# Check model with slack variables

# Stochastic Boot
println("########################")
slack_variables_boot = []
slack_variables_boot_index = []
for i in 1:sc
    foldername = "Stochastic_Bootstrap1_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        push!(slack_variables_boot_index,i)
        # Get placement from stochastic problem
        temp = get_solution_stochastic(foldername,"Stochastic_Solution",HPC_folder)
        cs = temp.cs
        model = second_stage_model_slack(cs, Deterministic_problem)
        set_silent(model) # removes terminal output
        set_time_limit_sec(model, 60 * 5) # 5 minutes to start with
        optimize!(model)
        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("slack model status: ", status)
        else
            push!(slack_variables_boot,[value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax])
            ])
        end
    end
end
# Slack variables
println("##########################")
println("Model indecies infeasible for Stochastic Boot: ", slack_variables_boot_index)
println("Slack variables for Stochastic Boot")
for i in 1:length(slack_variables_boot)
    for j in 1:length(slack_variables_boot[i])
        if any(slack_variables_boot[i][j] .!= 0.0)
            if j == 1
                println("In model $(slack_variables_boot_index[i]), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_boot[i][j])
            else
                println("Model $(slack_variables_boot_index[i]) and slack variable index $(j) are not zero")
                println(slack_variables_boot[i][j])
            end
        end
    end
end

# EVP Boot
#=
slack_variables_EVP_boot = []
slack_variables_EVP_boot_index = []
for i in 1:sc
    foldername = "EVP_Bootstrap1_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        push!(slack_variables_EVP_boot_index,i)
        # Get placement from stochastic problem
        temp = get_solution_deterministic(foldername,"EVP_Solution",HPC_folder)
        cs = temp.cs
        model = second_stage_model_slack(cs, Deterministic_problem)
        set_silent(model) # removes terminal output
        set_time_limit_sec(model, 60 * 5) # 5 minutes to start with
        optimize!(model)
        status = termination_status(model)
        if status != MOI.OPTIMAL
            println("slack model status: ", status)
        else
            push!(slack_variables_EVP_boot,[value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax])
            ])
        end
    end
end
# Slack variables
println("##########################")
println("Model indecies infeasible for EVP Boot: ", slack_variables_EVP_boot_index)
println("Slack variables for EVP Boot")
for i in 1:length(slack_variables_EVP_boot)
    for j in 1:length(slack_variables_EVP_boot[i])
        if any(slack_variables_EVP_boot[i][j] .!= 0.0)
            if j == 1
                println("In model $(slack_variables_EVP_boot_index[i]), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_EVP_boot[i][j])
            else
                println("Model $(slack_variables_EVP_boot_index[i]) and slack variable index $(j) are not zero")
                println(slack_variables_EVP_boot[i][j])
            end
        end
    end
end
=#
# Stochastic Gen 
slack_variables_gen = []
slack_variables_gen_index = []
for i in 1:sc
    foldername = "Stochastic_random_sampling_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    if (temp.status == "INFEASIBLE") # solve slack model
        push!(slack_variables_gen_index,i)
        # Get placement from stochastic problem
        temp = get_solution_stochastic(foldername,"Stochastic_Solution",HPC_folder)
        cs = temp.cs
        model = second_stage_model_slack(cs, Deterministic_problem)
        set_silent(model) # removes terminal output
        set_time_limit_sec(model, 60 * 5) # 5 minutes to start with
        optimize!(model)
        status = termination_status(model)
        println(typeof(status))
        if status != MOI.OPTIMAL
            println("slack model status: ", status)
        else
            push!(slack_variables_gen,[value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax])
            ])
        end
    end
end
# Slack variables
println("##########################")
println("Model indecies infeasible for Stochastic Gen: ", slack_variables_gen_index)
println("Slack variables for Stochastic Gen")
for i in 1:length(slack_variables_gen)
    for j in 1:length(slack_variables_gen[i])
        if any(slack_variables_gen[i][j] .!= 0.0)
            if j == 1
                println("In model $(slack_variables_gen_index[i]), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_gen[i][j])
            else
                println("Model $(slack_variables_gen_index[i]) and slack variable index $(j) are not zero")
                println(slack_variables_gen[i][j])
            end
        end
    end
end
# EVP Gen 
#=
slack_variables_EVP_gen = []
slack_variables_EVP_gen_index = []
for i in 1:sc
    foldername = "EVP_random_sampling_rep$(1)_sc$(scenarios[i])_unknown$(n_unknown[1])_time$(time_limit)"
    temp = get_solution_deterministic(foldername,
    "Fitted_Solution",HPC_folder)
    #println(temp.status)
    if (temp.status == "INFEASIBLE") # solve slack model
        push!(slack_variables_EVP_gen_index,i)
        # Get placement from stochastic problem
        temp = get_solution_deterministic(foldername,"EVP_Solution",HPC_folder)
        cs = temp.cs
        model = second_stage_model_slack(cs, Deterministic_problem)
        set_silent(model) # removes terminal output
        set_time_limit_sec(model, 60 * 5) # 5 minutes to start with
        optimize!(model)
        status = termination_status(model)
        println(typeof(status))
        if status != MOI.OPTIMAL
            println("slack model status: ", status)
            error("Model infeasible - do something to slack model")
        else
            push!(slack_variables_EVP_gen,[value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax])
            ])
        end
    end
end
# Slack variables
println("##########################")
println("Model indecies infeasible for EVP Gen: ", slack_variables_EVP_gen_index)
println("Slack variables for EVP Gen")
for i in 1:length(slack_variables_EVP_gen)
    for j in 1:length(slack_variables_EVP_gen[i])
        if any(abs(slack_variables_EVP_gen[i][j]) .!= 0.0)
            if j == 1
                println("In model $(slack_variables_EVP_gen_index[i]), the deck slack variable is not zero. Deck limit is violated")
                println(slack_variables_EVP_gen[i][j])
            else
                println("Model $(slack_variables_EVP_gen_index[i]) and slack variable index $(j) are not zero")
                println(slack_variables_EVP_gen[i][j])
            end
        end
    end
end
=#


# Looking into limiting factor - Deck 1
println("##########################")
n_decks = length(Deterministic_problem.vessel.decks)
deck_limits = [Deterministic_problem.vessel.decks[i].weight_limit for i in 1:n_decks]
println("Deterministic Solution weight on decks:")
#secu_var = [collect(filter(x -> x.QuantileNumber == i, df_secu_quan).Variance) for i in 1:n_quantiles]
Det_sol_PlacementDecks = [filter(x -> x.deck == i, Deterministic_Solution.cargo) for i in 1:n_decks]
for i in 1:n_decks
    println("Weight Limit on deck $(i): ", deck_limits[i])
    println("Weight on deck $(i): ", sum(Det_sol_PlacementDecks[i][j].weight for j in 1:length(Det_sol_PlacementDecks[i])))
end

xtemp = [[],[],[],[]]
deck_diff = [[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]]
for i in 1:repetitions
    for j in 1:sc
        for l in 1:n
            #if EVP_gen_fitted[i,j,l].status == "TIME_LIMIT" || EVP_gen_fitted[i,j,l].status == "OPTIMAL"
            #    EVP_gen_fitted_PlacementDecks = [filter(x -> x.deck == k, EVP_gen_fitted[i,j,l].cargo) for k in 1:n_decks]
            #    println("#################")
            #    push!(xtemp[1],scenarios[j])
            #    for k in 1:n_decks
            #        push!(deck_diff[k][1], deck_limits[k]-sum(EVP_gen_fitted_PlacementDecks[k][m].weight for m in 1:length(EVP_gen_fitted_PlacementDecks[k])))
            #        println("EVP gen fitted model $(i), scenarios $(scenarios[j]), unknown cargo $(n_unknown[l]) weight-difference on deck $(k): ", deck_limits[k]-sum(EVP_gen_fitted_PlacementDecks[k][m].weight for m in 1:length(EVP_gen_fitted_PlacementDecks[k])))
            #    end
            #else # not feasible, probably due to deck 1
            #    push!(xtemp[1],scenarios[j])
            #    EVP_gen_PlacementDecks = [filter(x -> x.deck == k,EVP_gen[i,j,l].cargo) for k in 1:n_decks]
                # calculate how much over limit - same as slack variables?
            #    for k in 1:n_decks
            #        actual_weight = [filter(x -> x.id == EVP_gen_PlacementDecks[k][i].id, Deterministic_problem.cargo.items)[1].weight for i in 1:length(EVP_gen_PlacementDecks[k])]
            #        push!(deck_diff[k][1], deck_limits[k]-sum(actual_weight))
            #    end
            #end
            #=
            if EVP_boot_fitted[i,j,l].status == "TIME_LIMIT" || EVP_boot_fitted[i,j,l].status == "OPTIMAL"
                EVP_boot_fitted_PlacementDecks = [filter(x -> x.deck == k, EVP_boot_fitted[i,j,l].cargo) for k in 1:n_decks]
                println("#################")
                push!(xtemp[2],scenarios[j])
                for k in 1:n_decks
                    push!(deck_diff[k][2], deck_limits[k]-sum(EVP_boot_fitted_PlacementDecks[k][m].weight for m in 1:length(EVP_boot_fitted_PlacementDecks[k])))
                    println("EVP boot fitted model $(i), scenarios $(scenarios[j]), unknown cargo $(n_unknown[l]) weight on deck $(k): ", deck_limits[k]-sum(EVP_boot_fitted_PlacementDecks[k][m].weight for m in 1:length(EVP_boot_fitted_PlacementDecks[k])))
                end
            else # not feasible, probably due to deck 1
                push!(xtemp[2],scenarios[j])
                EVP_boot_PlacementDecks = [filter(x -> x.deck == k,EVP_boot[i,j,l].cargo) for k in 1:n_decks]
                # calculate how much over limit - same as slack variables?
                for k in 1:n_decks
                    actual_weight = [filter(x -> x.id == EVP_boot_PlacementDecks[k][i].id, Deterministic_problem.cargo.items)[1].weight for i in 1:length(EVP_boot_PlacementDecks[k])]
                    push!(deck_diff[k][2], deck_limits[k]-sum(actual_weight))
                end
            end
            =#
            if Stochastic_gen_fitted[i,j,l].status == "TIME_LIMIT" || Stochastic_gen_fitted[i,j,l].status == "OPTIMAL"
                Stochastic_gen_fitted_PlacementDecks = [filter(x -> x.deck == k, Stochastic_gen_fitted[i,j,l].cargo) for k in 1:n_decks]
                println("#################")
                push!(xtemp[3],scenarios[j])
                for k in 1:n_decks
                    push!(deck_diff[k][3], deck_limits[k]-sum(Stochastic_gen_fitted_PlacementDecks[k][m].weight for m in 1:length(Stochastic_gen_fitted_PlacementDecks[k])))
                    println("Stochastic gen fitted model $(i), scenarios $(scenarios[j]), unknown cargo $(n_unknown[l]) weight on deck $(k): ", deck_limits[k]-sum(Stochastic_gen_fitted_PlacementDecks[k][m].weight for m in 1:length(Stochastic_gen_fitted_PlacementDecks[k])))
                end
            else # not feasible, probably due to deck 1
                push!(xtemp[3],scenarios[j])
                Stochastic_gen_PlacementDecks = [filter(x -> x.deck == k,Stochastic_gen[i,j,l].cargo) for k in 1:n_decks]
                # calculate how much over limit - same as slack variables?
                for k in 1:n_decks
                    actual_weight = [filter(x -> x.id == Stochastic_gen_PlacementDecks[k][i].id, Deterministic_problem.cargo.items)[1].weight for i in 1:length(Stochastic_gen_PlacementDecks[k])]
                    push!(deck_diff[k][3], deck_limits[k]-sum(actual_weight))
                end
            end
            if Stochastic_boot_fitted[i,j,l].status == "TIME_LIMIT" || Stochastic_boot_fitted[i,j,l].status == "OPTIMAL"
                Stochastic_boot_fitted_PlacementDecks = [filter(x -> x.deck == k, Stochastic_boot_fitted[i,j,l].cargo) for k in 1:n_decks]
                println("#################")
                push!(xtemp[4],scenarios[j])
                for k in 1:n_decks
                    push!(deck_diff[k][4], deck_limits[k]-sum(Stochastic_boot_fitted_PlacementDecks[k][m].weight for m in 1:length(Stochastic_boot_fitted_PlacementDecks[k])))
                    println("Stochastic boot fitted model $(i), scenarios $(scenarios[j]), unknown cargo $(n_unknown[l]) weight on deck $(k): ", deck_limits[k]-sum(Stochastic_boot_fitted_PlacementDecks[k][m].weight for m in 1:length(Stochastic_boot_fitted_PlacementDecks[k])))
                end
            else # not feasible, probably due to deck 1
                push!(xtemp[4],scenarios[j])
                Stochastic_boot_PlacementDecks = [filter(x -> x.deck == k,Stochastic_boot[i,j,l].cargo) for k in 1:n_decks]
                # calculate how much over limit - same as slack variables?
                for k in 1:n_decks
                    actual_weight = [filter(x -> x.id == Stochastic_boot_PlacementDecks[k][i].id, Deterministic_problem.cargo.items)[1].weight for i in 1:length(Stochastic_boot_PlacementDecks[k])]
                    push!(deck_diff[k][4], deck_limits[k]-sum(actual_weight))
                end
            end
        end
    end
end
# plot for decks
p = []
for i in 1:n_decks
    pl = plot(xtemp[3],deck_diff[i][3], label = "Sto-gen", marker = :circle, markersize = 4,
    xlabel="Scenarios", ylabel="Weight diff. (t)",
    title = "Unused Weight Capacity on Deck $(i)\n for different models")
    #plot!(xtemp[2],deck_diff[i][2], label = "EVP-boot", marker = :circle, markersize = 4)
    #plot!(xtemp[3],deck_diff[i][3], label = "Sto-gen", marker = :circle, markersize = 4)
    plot!(xtemp[4],deck_diff[i][4], label = "Sto-boot", marker = :circle, markersize = 4)
    savefig(plot_folder*"Deck$(i)_WeightDiff.png")
    push!(p,pl)
end
for i in p
    display(i)
end

