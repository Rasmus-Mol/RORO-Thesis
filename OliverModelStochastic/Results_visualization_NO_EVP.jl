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
                    filename,HPC_folders[l],problemname1,test_instance,problemname3)
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
                    filename,HPC_folders[l],problemname1,test_instance,problemname3)
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
# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end
gaps_boot = zeros(repetitions, sc,length(HPC_folders))
gaps_gen = zeros(repetitions, sc,length(HPC_folders))
cargo_loaded_boot = zeros(repetitions, sc,length(HPC_folders))
cargo_loaded_gen = zeros(repetitions, sc,length(HPC_folders))
ballast_water_boot_fitted = zeros(repetitions, sc,length(HPC_folders))
ballast_water_gen_fitted = zeros(repetitions, sc,length(HPC_folders))
ballast_water_boot = Array{Any}(undef,repetitions, sc,length(HPC_folders))
ballast_water_gen = Array{Any}(undef,repetitions, sc,length(HPC_folders))
objective_val_gen = zeros(repetitions, sc,length(HPC_folders))
objective_val_boot = zeros(repetitions, sc,length(HPC_folders))
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
            ballast_water_boot_fitted[i, j,l] = Stochastic_boot_fitted[i, j, 1,l].ballast_weight
            ballast_water_gen_fitted[i, j,l] = Stochastic_gen_fitted[i, j, 1,l].ballast_weight
            ballast_water_boot[i, j,l] = mean(getfield.(Stochastic_boot[i, j, 1,l].forces, :ballast_weight))
            ballast_water_gen[i, j,l] = mean(getfield.(Stochastic_gen[i, j, 1,l].forces, :ballast_weight))
            
            # 
            objective_val_gen[i, j,l] = Stochastic_gen[i, j, 1,l].objective
            objective_val_boot[i, j,l] = Stochastic_boot[i, j, 1,l].objective
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

temp = 0
for i in 1:repetitions
    temp += mean(Stochastic_problem_boot[i,1,1,4].cargo.items.total_weight)
end
temp/ repetitions

Stochastic_boot[9,5,1,8].cargo

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

# Objective value stuff

for j in 1:length(HPC_folders)
    reps = [objective_val_boot[:, i, j] for i in 1:length(scenarios)]
    #reps = [objective_val_gen[:, i, j] for i in 1:length(scenarios)]
    # 2) Compute mean, standard error, and 95% CI half-width for each N
    df = DataFrame(N=Int[], mean=Float64[], sem=Float64[], ci95=Float64[])
    for (i, N) in enumerate(scenarios)
        v = reps[i]
        μ = mean(v)
        σ = std(v, corrected=true)
        sem = σ / sqrt(length(v))
        hw = 1.96 * sem               # 95% CI half-width
        push!(df, (N, μ, sem, hw))
    end
    println("Confidence interval, instance $(j): ", df.ci95)
    # Plot mean ± 95% CI vs. N
    p = plot(df.N, df.mean;
        yerror=df.ci95,
        marker=:circle,
        xlabel="Number of scenarios",
        ylabel="Ojective value",
        title  = "Instance $(j)",
        legend=false,yformatter = :plain)
    display(p)
    savefig(p,plot_folder*"ObjectiveValue_vs_Scenarios_Instance_$(j)_boot.png")
    #savefig(p,plot_folder*"ObjectiveValue_vs_Scenarios_Instance_$(j)_gen.png")
end
# ballast water - fitted
for j in 1:length(HPC_folders)
    reps = [ballast_water_boot_fitted[:, i, j] for i in 1:length(scenarios)]
    #reps = [objective_val_gen[:, i, j] for i in 1:length(scenarios)]
    # 2) Compute mean, standard error, and 95% CI half-width for each N
    df = DataFrame(N=Int[], mean=Float64[], sem=Float64[], ci95=Float64[])
    for (i, N) in enumerate(scenarios)
        v = reps[i]
        μ = mean(v)
        σ = std(v, corrected=true)
        sem = σ / sqrt(length(v))
        hw = 1.96 * sem               # 95% CI half-width
        push!(df, (N, μ, sem, hw))
    end
    println("Confidence interval, instance $(j): ", df.ci95)
    # Plot mean ± 95% CI vs. N
    p = plot(df.N, df.mean;
        yerror=df.ci95,
        marker=:circle,
        xlabel="Number of scenarios",
        ylabel="Ballast water",
        title  = "Instance $(j)",
        legend=false,yformatter = :plain)
    display(p)
    #savefig(p,plot_folder*"ObjectiveValue_vs_Scenarios_Instance_$(j)_boot.png")
    #savefig(p,plot_folder*"ObjectiveValue_vs_Scenarios_Instance_$(j)_gen.png")
end
# ballast water 
for j in 1:length(HPC_folders)
    reps = [ballast_water_boot[:, i, j] for i in 1:length(scenarios)]
    #reps = [objective_val_gen[:, i, j] for i in 1:length(scenarios)]
    # 2) Compute mean, standard error, and 95% CI half-width for each N
    df = DataFrame(N=Int[], mean=Float64[], sem=Float64[], ci95=Float64[])
    for (i, N) in enumerate(scenarios)
        v = reps[i]
        μ = mean(v)
        σ = std(v, corrected=true)
        sem = σ / sqrt(length(v))
        hw = 1.96 * sem               # 95% CI half-width
        push!(df, (N, μ, sem, hw))
    end
    println("Confidence interval, instance $(j): ", df.ci95)
    # Plot mean ± 95% CI vs. N
    p = plot(df.N, df.mean;
        yerror=df.ci95,
        marker=:circle,
        xlabel="Number of scenarios",
        ylabel="Ballast water",
        title  = "Instance $(j)",
        legend=false,yformatter = :plain)
    display(p)
    #savefig(p,plot_folder*"ObjectiveValue_vs_Scenarios_Instance_$(j)_boot.png")
    #savefig(p,plot_folder*"ObjectiveValue_vs_Scenarios_Instance_$(j)_gen.png")
end


std(objective_val_boot, dims = 1)
p1 = plot(mean(std(objective_val_boot, dims = 1), dims = 3)[1,:,1], label = "",
ylabel = "Standard deviation", xlabel = "Number of scenarios", yformatter = :plain)
savefig(p1, plot_folder*"MeanStandardDeviation_vs_Scenarios_Boot.png")
p2 = plot(mean(std(objective_val_gen, dims = 1), dims = 3)[1,:,1], label = "",
ylabel = "Standard deviation", xlabel = "Number of scenarios", yformatter = :plain)
savefig(p2, plot_folder*"MeanStandardDeviation_vs_Scenarios_Gen.png")
# checking not all instances
idx_temp = vcat(2,4,)
temp_gen = objective_val_gen[:,:,idx_temp]
p1 = plot(mean(std(temp_gen, dims = 1), dims = 3)[1,:,1], label = "",
ylabel = "Standard deviation", xlabel = "Number of scenarios", yformatter = :plain)
savefig(p1, plot_folder*"MeanStandardDeviation_vs_Scenarios_Gen_instance_1_2.png")
temp_boot = objective_val_boot[:,:,idx_temp]
p2 = plot(mean(std(temp_boot, dims = 1), dims = 3)[1,:,1], label = "",
ylabel = "Standard deviation", xlabel = "Number of scenarios", yformatter = :plain)
savefig(p2, plot_folder*"MeanStandardDeviation_vs_Scenarios_Boot_instance_1_2.png")


Stochastic_problem_boot[1,1,1,1].cargo.items[1]
objective_val_boot[1,1,1]
Stochastic_boot[1,1,1,1].gap
avg_obj_boot = mean(objective_val_boot, dims = 1)
plot(avg_obj_boot[1,:,1])
size(avg_obj_boot)
getfield.(Stochastic_boot[:,1,1,1],:gap)
getfield.(Stochastic_boot[:,3,1,1],:gap)

# check if cargo is too heavy
for i in 1:length(HPC_folders)
    for j in 1:repetitions
        for k in 1:sc
            for l in 1:scenarios[k]
                if maximum([c.weight for c in filter(x->x.cargo_type_id == 1, Stochastic_problem_boot[j,k,1,i].cargo.items[l])]) > 40
                    println(
                        maximum([c.weight for c in filter(x->x.cargo_type_id == 1, Stochastic_problem_boot[j,k,1,i].cargo.items[l])])
                    )
                end
                if maximum([c.weight for c in filter(x->x.cargo_type_id == 4, Stochastic_problem_boot[j,k,1,i].cargo.items[l])]) > 99.5
                    println(
                        maximum([c.weight for c in filter(x->x.cargo_type_id == 4, Stochastic_problem_boot[j,k,1,i].cargo.items[l])])
                    )
                end
            end
        end
    end
end

# Scenarios
plot_folder_scenarios = "Plots/Data/Scenarios/"

cur_path = @__DIR__
file_path = joinpath(cur_path,"data","CargoWeights.csv")
#joinpath(cur_path,"..","..","data","CargoWeights.csv")
df = load_Weight_Variance_data(file_path)
df_secu, df_trailer = weight_difference_info(df, false)
sort_df_secu = sort(df_secu, :CountBookedWeight)

p1 = plot(sort_df_secu.CountBookedWeight, sort_df_secu.Variance,
title = "Weight difference for Secu-boxes",xlabel = "Booked weight (t.)", ylabel = "Weight difference (t.)",
yformatter = :plain, xformatter = :plain)
savefig(p1, "Plots/Data/BookedWeight_vs_WeightDiff_Secu_.png")
sort_df_trailer = sort(df_trailer, :CountBookedWeight)
p2 = plot(sort_df_trailer.CountBookedWeight, sort_df_trailer.Variance,
title = "Weight difference for Trucks",xlabel = "Booked weight (t.)", ylabel = "Weight difference (t.)",
yformatter = :plain, xformatter = :plain)
savefig(p2, "Plots/Data/BookedWeight_vs_WeightDiff_Trailer_.png")

n_quantiles = 4
df_secu, q_arr_secu = seperate_data_into_quantiles(df_secu, n_quantiles, false)
df_trailer, q_arr_trailer = seperate_data_into_quantiles(df_trailer, n_quantiles, false)
q_arr_secu = q_arr_secu ./ 1000 # convert to tons
q_arr_trailer = q_arr_trailer ./ 1000 # convert to tons
secu_var = [collect(filter(x -> x.QuantileNumber == i, df_secu).Variance) for i in 1:n_quantiles]
trailer_var = [collect(filter(x -> x.QuantileNumber == i, df_trailer).Variance) for i in 1:n_quantiles]


# Comparing instance 1 and 2 in SGM 1 and 2

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

n_quantiles = 4
df = load_Weight_Variance_data(file_path)
df_secu, df_trailer = weight_difference_info(df, false)
df_secu, q_arr_secu = seperate_data_into_quantiles(df_secu, n_quantiles, false)
df_trailer, q_arr_trailer = seperate_data_into_quantiles(df_trailer, n_quantiles, false)
q_arr_secu = q_arr_secu ./ 1000 # convert to tons
q_arr_trailer = q_arr_trailer ./ 1000 # convert to tons

secu_var = [collect(filter(x -> x.QuantileNumber == i, df_secu).Variance) for i in 1:n_quantiles]./1000
trailer_var = [collect(filter(x -> x.QuantileNumber == i, df_trailer).Variance) for i in 1:n_quantiles]./1000
for i in 1:length(secu_var)
    println("##############")
    println("mean secu: ", mean(secu_var[i]))
    println("Max secu: ", maximum(secu_var[i]))
    println("Min secu: ", minimum(secu_var[i]))
    println("mean Trucks: ", mean(trailer_var[i]))
    println("Max Trucks: ", maximum(trailer_var[i]))
    println("Min Trucks: ", minimum(trailer_var[i]))
end


using Distributions
using StatsPlots       # for histogram, pdf overlay, qqplot
using HypothesisTests  # for ks test
using StatsBase
using Statistics
using Optim

pvals = Vector{Vector{Float64}}(undef, length(secu_var))
for i in 1:length(secu_var)
    data_og = secu_var[i]
    # discrete
    lo, hi = minimum(data_og), maximum(data_og)
    lo = round(lo/0.5)*0.5
    hi = round(hi/0.5)*0.5
    edges = lo:0.5:hi   # bins [lo,lo+0.5), [lo+0.5,lo+1.0), …, [hi-0.5,hi]
    h = fit(Histogram, data_og, edges)
    counts = h.weights              # integer counts
    k = length(counts)
    # Discrete uniform over these k bins:
    probs = fill(1/k, k)
    χ2 = ChisqTest(counts, probs)
    println(χ2)
    pvals[1] = [pvalue(χ2)]

    # Student t-tests
    # 1) Compute location & scale from the raw data
    μ, σ = mean(data_og), std(data_og)
    # 2) Standardize
    z = (data_og .- μ) ./ σ
    # 3) Define the negative log-likelihood as a function of ν
    function negloglik(ν)
        ν <= 0 && return Inf         # ν must be positive
        return -sum(logpdf.(TDist(ν), z))
    end
    # 4) Optimize ν over a reasonable range, say [1e-3, 200]
    res = optimize(negloglik, 1e-3, 200.0)
    ν̂  = Optim.minimizer(res)
    # 5) Reconstruct the fitted LocationScale t-distribution
    fit_td = LocationScale(μ, σ, TDist(ν̂))
    # Now you can do goodness-of-fit tests on `fit_td`
    #println("Fitted t: μ=$(round(μ,digits=3)), σ=$(round(σ,digits=3)), ν=$(round(ν̂,digits=3))")
    # Recreate your fitted t‐distribution
    fit_td = LocationScale(μ, σ, TDist(ν̂ ))
    # Anderson–Darling test
    ad = OneSampleADTest(data_og, fit_td)
    # Print the test result
    println(ad)
    push!(pvals,pvalue(ad))
    pvalue(ad)
end


data = data_og .+ (rand(length(data)) .- 0.5).*1e-6

# 1. Fit candidate distributions
dists = Dict{String,UnivariateDistribution}()
dists["Normal"]     = fit(Normal, data)
dists["Exponential"]= fit(Exponential, data .- minimum(data) .+ eps()) 
dists["LogNormal"]  = fit(LogNormal, data[data .> 0])  # only positives

# 1) Run the Shapiro–Wilk test
sw = ShapiroWilkTest(data_og)

# 2) Extract the W statistic and p‐value
W  = round(sw.W, digits=4)
p  = round(pvalue(sw),    digits=4)
sw.W
println("Shapiro–Wilk Test:")
println("  W = $W")
println("  p = $p")

# 2. Histogram + PDF overlay
histogram(data; bins=30, normalize=true, alpha=0.4, label="data")
x = range(minimum(data), stop=maximum(data), length=200)
for (name, dist) in dists
    pdf_vals = pdf.(dist, x)
    plot!(x, pdf_vals; lw=2, label="$name fit")
end
xlabel!("Value"); ylabel!("Density")
title!("Data with Fitted PDFs")

# 3. Q–Q plots
for (name, dist) in dists
    display(qqplot(dist, data; title="Q–Q plot vs $name", legend=false))
end

# 4. KS test
println("Kolmogorov–Smirnov tests:")
for (name, dist) in dists
    ks = ApproximateOneSampleKSTest(data, dist)
    D = round(ks.δ, digits=3)          # overall KS statistic
    p = round(pvalue(ks), digits=3)     # p-value
    println("  $name: D=$D, p=$p")
end
ks = ApproximateOneSampleKSTest(data, dists["Normal"])
fieldnames(typeof(ks))
dump(ks)
statistic(ks)

# 5. AIC comparison
function aic(dist::UnivariateDistribution, data)
    loglik = sum(logpdf.(dist, data))
    k = length(fieldnames(typeof(dist)))  # approximate # of params
    return 2k - 2loglik
end

println("\nAIC scores:")
for (name, dist) in dists
    println("  $name: AIC=$(round(aic(dist, data), digits=1))")
end











