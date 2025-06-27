include("packages_and_files.jl")

HPC_folders = [
    "Finlandia_mixed_light_60_SlackDeck1_19_06_20",
    "Finlandia_mixed_light_100_SlackDeck1_19_06_20",
    "Finlandia_mixed_heavy_60_SlackDeck1_19_06_20",
    "Finlandia_mixed_heavy_100_SlackDeck1_19_06_21",
    "Finlandia_no_cars_light_60_SlackDeck1_19_06_21",
    "Finlandia_no_cars_light_100_SlackDeck1_19_06_21",
    "Finlandia_no_cars_heavy_60_SlackDeck1_19_06_21",
    "Finlandia_no_cars_heavy_100_SlackDeck1_19_06_21",
]
problemname1, problemname3 = "finlandia", "hazardous"
repetitions = 5
slack_fraction = [[1.1, 1, 1], [1.2, 1, 1], [1.3, 1, 1], [1.4, 1, 1], [1.5, 1, 1]]
scenarios = [10]

Deterministic_problem = Array{Any}(undef, length(HPC_folders))#, length(slack_fraction))
Deterministic_Solution = Array{Any}(undef, length(HPC_folders), length(slack_fraction))

Stochastic_gen = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
Stochastic_gen_fitted = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
Stochastic_boot = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
Stochastic_boot_fitted = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))

boot_fitted_slacked = Array{Any}(nothing, length(HPC_folders), repetitions, length(slack_fraction))
gen_fitted_slacked = Array{Any}(nothing, length(HPC_folders), repetitions, length(slack_fraction))
slack_boot = Array{Any}(nothing, length(HPC_folders), repetitions, length(slack_fraction))
slack_gen = Array{Any}(nothing, length(HPC_folders), repetitions, length(slack_fraction))

Stochastic_problem_gen = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
Stochastic_problem_boot = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))


# EVP solutions
EVP_gen = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
EVP_gen_pro = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
EVP_gen_matrix = Array{Any}(undef, length(HPC_folders), length(slack_fraction), 4)

EVP_boot = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
EVP_boot_pro = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
EVP_boot_matrix = Array{Any}(undef, length(HPC_folders), length(slack_fraction), 4)
# For later to get results
function get_obj_of_evp_2D(filename::String, foldername::String, HPC_folder::String)
    temp = open(joinpath("Results", HPC_folder, foldername, filename * ".json"), "r") do file
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

# load results
for i in 1:length(HPC_folders)
    test_instance = Finlandia_test[i]
    HPC_folder = HPC_folders[i]
    rep, sc, n_unknown, time_limit, note = get_HPC_data(HPC_folder)
    #n = length(n_unknown)
    Deterministic_problem[i] = load_data(problemname1, test_instance, problemname3)
    for j in 1:length(slack_fraction)
        Deterministic_Solution[i, j] = get_solution_deterministic("Finlandia_deterministic",
            "Deterministic_Solution_slack_$(slack_fraction[j])", HPC_folder)
        # EVP_matrix: obj,ballast, cargo, infeasibility
        foldername = "EVP_gen_info_slack_$(slack_fraction[j])"
        EVP_gen_matrix[i, j, 1] = get_obj_of_evp_2D("Objective_values_gen",
            foldername, HPC_folder)
        EVP_gen_matrix[i, j, 2] = get_obj_of_evp_2D("Ballast_used_gen",
            foldername, HPC_folder)
        EVP_gen_matrix[i, j, 3] = get_obj_of_evp_2D("Cargo_loaded_gen",
            foldername, HPC_folder)
        EVP_gen_matrix[i, j, 4] = get_obj_of_evp_2D("inf_index_gen",
            foldername, HPC_folder)
        # EVP boot
        foldername = "EVP_boot_info_slack_$(slack_fraction[j])"
        EVP_boot_matrix[i, j, 1] = get_obj_of_evp_2D("Objective_values_boot",
            foldername, HPC_folder)
        EVP_boot_matrix[i, j, 2] = get_obj_of_evp_2D("Ballast_used_boot",
            foldername, HPC_folder)
        EVP_boot_matrix[i, j, 3] = get_obj_of_evp_2D("Cargo_loaded_boot",
            foldername, HPC_folder)
        EVP_boot_matrix[i, j, 4] = get_obj_of_evp_2D("inf_index_boot",
            foldername, HPC_folder)
        for k in 1:repetitions
            # EVP Gen
            foldername = "EVP_random_sampling_rep$(k)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)_slack_$(slack_fraction[j])"
            filename = "EVP_Solution"
            EVP_gen[i, k, j] = get_solution_deterministic(foldername,
                filename, HPC_folder)
            filename = "EVP_Problem"
            EVP_gen_pro[i, k, j] = get_deterministic_problem(foldername,
                filename, HPC_folder, problemname1, test_instance, problemname3)

            # EVP Boot
            foldername = "EVP_Bootstrap1_rep$(k)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)_slack_$(slack_fraction[j])"
            filename = "EVP_Solution"
            EVP_boot[i, k, j] = get_solution_deterministic(foldername,
                filename, HPC_folder)
            filename = "EVP_Problem"
            EVP_boot_pro[i, k, j] = get_deterministic_problem(foldername,
                filename, HPC_folder, problemname1, test_instance, problemname3)

            # Stochastic Gen
            foldername = "Stochastic_random_sampling_rep$(k)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)_slack_$(slack_fraction[j])"
            filename = "Stochastic_Problem"
            # problem
            # Didn't save gen problem
            #Stochastic_problem_gen[i,k,j] = get_stochastic_problem(foldername,
            #    filename,HPC_folder,problemname1,test_instance,problemname3)
            # solution
            filename = "Stochastic_Solution"
            Stochastic_gen[i, k, j] = get_solution_stochastic(foldername,
                filename, HPC_folder)
            Stochastic_gen_fitted[i, k, j] = get_solution_deterministic(foldername,
                "Fitted_Solution", HPC_folder)
            # Load slacked
            if Stochastic_gen_fitted[i, k, j].status != "OPTIMAL"
                gen_fitted_slacked[i, k, j] = get_solution_deterministic(foldername,
                    "Fitted_Solution_slacked", HPC_folder)
                slack_gen[i, k, j] = get_slack(foldername, "Fitted_Solution", HPC_folder)
            end

            # Stochastic Boot
            foldername = "Stochastic_Bootstrap1_rep$(k)_sc$(scenarios[1])_unknown$(n_unknown[1])_time$(time_limit)_slack_$(slack_fraction[j])"
            filename = "Stochastic_Problem"
            # problem
            # i:folders, j:slack, k:rep
            Stochastic_problem_boot[i, k, j] = get_stochastic_problem(foldername,
                filename, HPC_folder, problemname1, test_instance, problemname3)
            # solution
            filename = "Stochastic_Solution"
            Stochastic_boot[i, k, j] = get_solution_stochastic(foldername,
                filename, HPC_folder)
            Stochastic_boot_fitted[i, k, j] = get_solution_deterministic(foldername,
                "Fitted_Solution", HPC_folder)
            # Load slacked
            if Stochastic_boot_fitted[i, k, j].status != "OPTIMAL"
                boot_fitted_slacked[i, k, j] = get_solution_deterministic(foldername,
                    "Fitted_Solution_slacked", HPC_folder)
                slack_boot[i, k, j] = get_slack(foldername, "Fitted_Solution", HPC_folder)
            end
        end
    end
end
# results for table
gaps_det = zeros(length(HPC_folders), length(slack_fraction))
cargo_laoded_det = zeros(length(HPC_folders), length(slack_fraction))
ballast_water_det = zeros(length(HPC_folders), length(slack_fraction))
gaps_boot = zeros(length(HPC_folders), repetitions, length(slack_fraction))
gaps_gen = zeros(length(HPC_folders), repetitions, length(slack_fraction))
cargo_loaded_boot = zeros(length(HPC_folders), repetitions, length(slack_fraction))
cargo_loaded_gen = zeros(length(HPC_folders), repetitions, length(slack_fraction))
ballast_water_boot_fitted = zeros(length(HPC_folders), repetitions, length(slack_fraction))
ballast_water_gen_fitted = zeros(length(HPC_folders), repetitions, length(slack_fraction))
ballast_water_boot = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
ballast_water_gen = Array{Any}(undef, length(HPC_folders), repetitions, length(slack_fraction))
objective_val_gen = zeros(length(HPC_folders), repetitions, length(slack_fraction))
objective_val_boot = zeros(length(HPC_folders), repetitions, length(slack_fraction))
for i in 1:length(HPC_folders)
    for j in 1:repetitions
        for k in 1:length(slack_fraction)
            if Stochastic_boot_fitted[i, j, k].status != "OPTIMAL"
                println("Boot - status: ", Stochastic_boot_fitted[i, j, k].status)
                println("test: $(i), rep: $(j), slack: $(k)")
            end
            if Stochastic_gen_fitted[i, j, k].status != "OPTIMAL"
                println("gen - status: ", Stochastic_gen_fitted[i, j, k].status)
                println("test: $(i), rep: $(j), slack: $(k)")
            end

            gaps_boot[i, j, k] = round(Stochastic_boot[i, j, k].gap * 100, digits=3)
            gaps_gen[i, j, k] = round(Stochastic_gen[i, j, k].gap * 100, digits=3)
            cargo_loaded_boot[i, j, k] = Int(Stochastic_boot[i, j, k].n_cargo_loaded)
            cargo_loaded_gen[i, j, k] = Int(Stochastic_gen[i, j, k].n_cargo_loaded)
            ballast_water_boot_fitted[i, j, k] = Stochastic_boot_fitted[i, j, k].ballast_weight
            ballast_water_gen_fitted[i, j, k] = Stochastic_gen_fitted[i, j, k].ballast_weight
            ballast_water_boot[i, j, k] = mean(getfield.(Stochastic_boot[i, j, k].forces, :ballast_weight))
            ballast_water_gen[i, j, k] = mean(getfield.(Stochastic_gen[i, j, k].forces, :ballast_weight))

            # 
            objective_val_gen[i, j, k] = Stochastic_gen[i, j, k].objective
            objective_val_boot[i, j, k] = Stochastic_boot[i, j, k].objective

            #
            gaps_det[i, k] = round(Deterministic_Solution[i, k].gap * 100, digits=3)
            cargo_laoded_det[i, k] = Deterministic_Solution[i, k].n_cargo_loaded
            ballast_water_det[i, k] = Deterministic_Solution[i, k].ballast_weight
        end
    end
end

for l in 1:length(slack_fraction)
    println("Slack: ", l)
    #pretty_table_1 = Array{Any}(undef,8, 5)
    pretty_table_1 = zeros(8,10)
    for i in 1:length(HPC_folders)
        #println("##############")
        #println("Instance: $(i)")
        #println("Det gap: ", gaps_det[i,:])
        #println("Cargo loaded: ",cargo_laoded_det[i,:])
        #println("Ballast water: ",ballast_water_det[i,:])
        #println("Det Status: ", Deterministic_Solution[i].status)
        #println(round(Deterministic_Solution[i].ballast_weight,digits=2))
        #println("Gap: ",Deterministic_Solution[i].gap)
        #println("Cargo loaded: ", Deterministic_Solution[i].n_cargo_loaded,"/",
        #Deterministic_Solution[i].n_cargo_total)
        #for j in 1:length(scenarios)

        #println("Gap: ", mean(gaps_gen[i,:,1], dims = 1))
        pretty_table_1[i, 1] = mean(gaps_gen[i, :, l])
        pretty_table_1[i, 2] = std(gaps_gen[i, :, l])
        pretty_table_1[i, 3] = mean(cargo_loaded_gen[i, :, l])
        pretty_table_1[i, 4] = mean(ballast_water_gen_fitted[i, :, l])
        #println("SD: ", std(gaps_gen[i,:,1], dims = 1))
        #println("Cargo loaded: ", mean(cargo_loaded_gen[i,:,1], dims = 1))
        #println("Ballast water: ", mean(ballast_water_gen[i,:,:],dims = 1))
        #println("Ballast water: ", mean(ballast_water_gen_fitted[i,:,1],dims = 1))
        temp = count(x -> x == nothing, gen_fitted_slacked[i, :, l])
        pretty_table_1[i, 5] = temp
        pretty_table_1[i, 6] = mean(gaps_boot[i, :, l])
        pretty_table_1[i, 7] = std(gaps_boot[i, :, l])
        pretty_table_1[i, 8] = mean(cargo_loaded_boot[i, :, l])
        pretty_table_1[i, 9] = mean(ballast_water_boot_fitted[i, :, l])
        #println("SD: ", std(gaps_gen[i,:,1], dims = 1))
        #println("Cargo loaded: ", mean(cargo_loaded_gen[i,:,1], dims = 1))
        #println("Ballast water: ", mean(ballast_water_gen[i,:,:],dims = 1))
        #println("Ballast water: ", mean(ballast_water_gen_fitted[i,:,1],dims = 1))
        temp = count(x -> x == nothing, boot_fitted_slacked[i, :, l])
        pretty_table_1[i, 10] = temp
        #header = [""
        #=
        println("Gap: ", mean(gaps_boot[i,:,:], dims = 1))
        println("SD: ", std(gaps_boot[i,:,:], dims = 1))
        println("Cargo loaded: ", mean(cargo_loaded_boot[i,:,:], dims = 1))
        #println("Ballast water: ", mean(ballast_water_boot[i,:,:],dims = 2))
        println("Ballast water: ", mean(ballast_water_boot_fitted[i,:,:],dims = 1))
        temp = zeros(length(slack_fraction))
        for j in 1:length(slack_fraction)
            temp[j] = count(x-> x== nothing, boot_fitted_slacked[i,:,j])
        end
        println("Feasible, sc: ", temp)
        =#
        #end
    end
    header = ["Col1", "Col2", "Col3", "Col4", "Col5"]
    pretty_table(pretty_table_1)
    #println(pretty_table)
end


# Comparing ballast water used to deterministic
ballast_w = Array{Any}(undef, length(HPC_folders), length(slack_fraction))
for i in 1:length(slack_fraction)
    for j in 1:length(HPC_folders)
        temp = getfield.(Stochastic_boot_fitted[j, :, i], :ballast_weight)
        temp = filter(x -> x != 0.0, temp) # remove zero ballast water
        if length(temp) > 0
            ballast_w[j, i] = mean(temp) - Deterministic_Solution[j, i].ballast_weight
        else
            ballast_w[j, i] = [0.0]
        end
    end
    #println(round.(filter(z -> z != [0.0], ballast_w[:, i]), digits=3))
    x = [k for k in 1:length(HPC_folders) if ballast_w[k, i] != [0.0]]
    #println(x)
    p = plot(x, round.(filter(z -> z != [0.0], ballast_w[:, i]), digits=3), label="",
        title="Mean ballast water difference between \n
      stochastic and deterministic solution",
        xlabel="Instance", ylabel="Difference in ballast water")
    display(p)
    savefig(p, plot_folder * "/Ballast_water_difference_slack_$(i).png")
end

# VVS & EVPI solutions
function get_obj_val(cargo_loaded::Float64, det_pro::StowageProblem, ballast_used::Float64)
    vessel = det_pro.vessel
    n_ballast_tanks = length(vessel.ballast_tanks)
    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
    #haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
    #cost = [CSC for c in cargo]
    #cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
    #area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]
    return ballast_used - CSC * cargo_loaded
end
corrected_obj_gen = Array{Any}(nothing, length(HPC_folders), repetitions, length(slack_fraction))
corrected_obj_boot = Array{Any}(nothing, length(HPC_folders), repetitions, length(slack_fraction))
EEV_gen = zeros(length(HPC_folders), repetitions, length(slack_fraction))
EEV_boot = zeros(length(HPC_folders), repetitions, length(slack_fraction))
sto_obj_gen = zeros(length(HPC_folders), repetitions, length(slack_fraction))
sto_obj_boot = zeros(length(HPC_folders), repetitions, length(slack_fraction))
sto_obj_fitted_gen = zeros(length(HPC_folders), repetitions, length(slack_fraction))
sto_obj_fitted_boot = zeros(length(HPC_folders), repetitions, length(slack_fraction))
VSS_gen = zeros(length(HPC_folders), repetitions, length(slack_fraction))
VSS_boot = zeros(length(HPC_folders), repetitions, length(slack_fraction))
EVPI_gen = zeros(length(HPC_folders), repetitions, length(slack_fraction))
EVPI_boot = zeros(length(HPC_folders), repetitions, length(slack_fraction))
# EVP_matrix: obj, ballast, cargo, infeasibility
for i in 1:length(HPC_folders)
    for j in 1:length(slack_fraction)
        for k in 1:repetitions
            # fix objective value
            temp_gen = zeros(10)
            temp_boot = zeros(10)
            for l in 1:10
                temp_gen[l] = get_obj_val(Float64(EVP_gen_matrix[i,j,3][k,l]), 
                Deterministic_problem[i], Float64(EVP_gen_matrix[i,j,2][k,l]))
                temp_boot[l] = get_obj_val(Float64(EVP_boot_matrix[i,j,3][k,l]), 
                Deterministic_problem[i], Float64(EVP_boot_matrix[i,j,2][k,l]))
            end
            corrected_obj_gen[i,k,j] = temp_gen
            corrected_obj_boot[i,k,j] = temp_boot
            # Finding EEV
            EEV_gen[i,k,j] = mean(temp_gen)
            EEV_boot[i,k,j] = mean(temp_boot)
            # Saving stochastic solution
            sto_obj_gen[i,k,j] = Stochastic_gen[i,k,j].objective
            sto_obj_boot[i,k,j] = Stochastic_boot[i,k,j].objective
            if Stochastic_gen_fitted[i,k,j].objective == Inf # infeasible
                sto_obj_fitted_gen[i,k,j] = get_obj_val(
                    Float64(Stochastic_gen_fitted[i,k,j].n_cargo_loaded),
                    Deterministic_problem[i],
                    Float64(Stochastic_gen_fitted[i,k,j].ballast_weight)
                )
            else # feasible
                sto_obj_fitted_gen[i,k,j] = Stochastic_gen_fitted[i,k,j].objective
            end
            if Stochastic_boot_fitted[i,k,j].objective == Inf # infeasible
                sto_obj_fitted_boot[i,k,j] = get_obj_val(
                    Float64(Stochastic_boot_fitted[i,k,j].n_cargo_loaded),
                    Deterministic_problem[i],
                    Float64(Stochastic_boot_fitted[i,k,j].ballast_weight)
                )
            else # feasible
                sto_obj_fitted_boot[i,k,j]= Stochastic_boot_fitted[i,k,j].objective
            end
            # VSS - using original stochastic solution
            VSS_gen[i,k,j] = EEV_gen[i,k,j] - sto_obj_gen[i,k,j]
            VSS_boot[i,k,j] = EEV_boot[i,k,j] - sto_obj_boot[i,k,j]
            # EVPI - using fitted stochastic solution
            EVPI_gen[i,k,j] = sto_obj_fitted_gen[i,k,j] - Deterministic_Solution[i,j].objective
            EVPI_boot[i,k,j] = sto_obj_fitted_boot[i,k,j] - Deterministic_Solution[i,j].objective
        end
    end
end
EVPI_gen[6,:,:]
Stochastic_boot_fitted[6,k,j]
EVPI_gen[:,:,4]
EVPI_boot[:,:,4]
getfield.(Stochastic_boot_fitted[1,:,:],:ballast_weight)


plot_folder = "Plots/Results/SlackDeck/"
# VSS
p1 = boxplot([VSS_gen[1,:,4], VSS_gen[3,:,4],VSS_gen[5,:,4],
    #VSS_boot[6,:,i], 
    VSS_gen[7,:,4],VSS_gen[8,:,4]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:5, ["1", "3", "5", "7", "8"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "VSS",           # axis labels
    title = "Value of Stochastic Solution - SGM 1")
savefig(p1,plot_folder*"/VSS_gen_slack_4_1.png")
p2 = boxplot([VSS_gen[2,:,4],VSS_gen[2,:,4]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:2, ["2", "4"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "VSS",           # axis labels
    title = "Value of Stochastic Solution - SGM 1")
savefig(p2,plot_folder*"/VSS_gen_slack_4_2.png")
p3 = boxplot([VSS_gen[6,:,4]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:1, ["6"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "VSS",           # axis labels
    title = "Value of Stochastic Solution - SGM 1")
savefig(p3,plot_folder*"/VSS_gen_slack_4_3.png")
p4 = boxplot([VSS_gen[6,1:4,4]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:1, ["6"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "VSS",           # axis labels
    title = "Value of Stochastic Solution - SGM 1")
savefig(p4,plot_folder*"/VSS_gen_slack_4_3_2.png")
# EVPI
# Gen
EVPI_gen[6,:,4]
p1 = boxplot([EVPI_gen[1,:,4], EVPI_gen[3,:,4],EVPI_gen[5,:,4],
#EVPI_gen[6,:,i], 
EVPI_gen[7,:,4],EVPI_gen[8,:,4]]; 
#labels = ["Instance 1", "Instance 3", "Instance 5 C",
#"Instance 6", "Instance 7"],    # custom labels for each box
labels = "",
xticks = (1:5, ["1", "3", "5", "7", "8"]),
#color = [:teal, :orange, :purple],            # fill colors for the boxes
xlabel = "Instance", 
ylabel = "EVPI",           # axis labels
title = "Expected Value of Perfect Information - SGM 1")
savefig(p1,plot_folder*"/EVPI_gen_slack_4_1.png")
p2 = boxplot([EVPI_gen[2,:,4],EVPI_gen[2,:,4]]; 
#labels = ["Instance 1", "Instance 3", "Instance 5 C",
#"Instance 6", "Instance 7"],    # custom labels for each box
labels = "",
xticks = (1:2, ["2", "4"]),
#color = [:teal, :orange, :purple],            # fill colors for the boxes
xlabel = "Instance", 
ylabel = "EVPI",           # axis labels
title = "Expected Value of Perfect Information - SGM 1")
savefig(p2,plot_folder*"/EVPI_gen_slack_4_2.png")
p3 = boxplot([EVPI_gen[6,:,4]]; 
#labels = ["Instance 1", "Instance 3", "Instance 5 C",
#"Instance 6", "Instance 7"],    # custom labels for each box
labels = "",
xticks = (1:1, ["6"]),
#color = [:teal, :orange, :purple],            # fill colors for the boxes
xlabel = "Instance", 
ylabel = "EVPI",           # axis labels
title = "Expected Value of Perfect Information - SGM 1")
savefig(p3,plot_folder*"/EVPI_gen_slack_4_3.png")
p4 = boxplot([EVPI_gen[6,1:4,4]]; 
#labels = ["Instance 1", "Instance 3", "Instance 5 C",
#"Instance 6", "Instance 7"],    # custom labels for each box
labels = "",
xticks = (1:1, ["6"]),
#color = [:teal, :orange, :purple],            # fill colors for the boxes
xlabel = "Instance", 
ylabel = "EVPI",           # axis labels
title = "Expected Value of Perfect Information - SGM 1")
savefig(p4,plot_folder*"/EVPI_gen_slack_4_3_2.png")

# VSS
for i in 1:length(slack_fraction)
    p1 = boxplot([VSS_boot[1,:,i], VSS_boot[3,:,i],VSS_boot[5,:,i],
    VSS_boot[6,:,i], VSS_boot[7,:,i],VSS_boot[8,:,i]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:6, ["1", "3", "5", "6", "7", "8"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "VSS",           # axis labels
    title = "Value of Stochastic Solution - SGM 2")
    savefig(p1,plot_folder*"/VSS_boot_slack_$(i)_1.png")
    p1 = boxplot([VSS_boot[2,:,i],VSS_boot[2,:,i]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:2, ["2", "4"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "VSS",           # axis labels
    title = "Value of Stochastic Solution - SGM 2")
    savefig(p1,plot_folder*"/VSS_boot_slack_$(i)_2.png")
    #
    p1 = boxplot([VSS_gen[1,:,i], VSS_gen[3,:,i],VSS_gen[5,:,i],
    VSS_gen[6,:,i], VSS_gen[7,:,i],VSS_gen[8,:,i]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:6, ["1", "3", "5", "6", "7", "8"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "VSS",           # axis labels
    title = "Value of Stochastic Solution - SGM 1")
    savefig(p1,plot_folder*"/VSS_gen_slack_$(i)_1.png")
    p1 = boxplot([VSS_gen[2,:,i],VSS_gen[2,:,i]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:2, ["2", "4"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "VSS",           # axis labels
    title = "Value of Stochastic Solution - SGM 1")
    savefig(p1,plot_folder*"/VSS_gen_slack_$(i)_2.png")
end
# EVPI
for i in 1:length(slack_fraction)
    p1 = boxplot([EVPI_boot[1,:,i], EVPI_boot[3,:,i],EVPI_boot[5,:,i],
    EVPI_boot[6,:,i], EVPI_boot[7,:,i],EVPI_boot[8,:,i]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:6, ["1", "3", "5", "6", "7", "8"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "EVPI",           # axis labels
    title = "Expected Value of Perfect Information - SGM 2")
    savefig(p1,plot_folder*"/EVPI_boot_slack_$(i)_1.png")
    p1 = boxplot([EVPI_boot[2,:,i],EVPI_boot[2,:,i]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:2, ["2", "4"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "EVPI",           # axis labels
    title = "Expected Value of Perfect Information - SGM 2")
    savefig(p1,plot_folder*"/EVPI_boot_slack_$(i)_2.png")
    # Gen
    p1 = boxplot([EVPI_gen[1,:,i], EVPI_gen[3,:,i],EVPI_gen[5,:,i],
    EVPI_gen[6,:,i], EVPI_gen[7,:,i],EVPI_gen[8,:,i]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:6, ["1", "3", "5", "6", "7", "8"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "EVPI",           # axis labels
    title = "Expected Value of Perfect Information - SGM 1")
    savefig(p1,plot_folder*"/EVPI_gen_slack_$(i)_1.png")
    p1 = boxplot([EVPI_gen[2,:,i],EVPI_gen[2,:,i]]; 
    #labels = ["Instance 1", "Instance 3", "Instance 5 C",
    #"Instance 6", "Instance 7"],    # custom labels for each box
    labels = "",
    xticks = (1:2, ["2", "4"]),
    #color = [:teal, :orange, :purple],            # fill colors for the boxes
    xlabel = "Instance", 
    ylabel = "EVPI",           # axis labels
    title = "Expected Value of Perfect Information - SGM 1")
    savefig(p1,plot_folder*"/EVPI_gen_slack_$(i)_2.png")
end




############### Ballast water comparison with origanl tests #######################
og_ball_gen = [592.32, 739.77, 102.33, 364.92, 481.63, 482.97,40.34]
og_ball_boot = [592.32,521.48, 93.11, 337.14,481.63, 342.89,46.79]
new_ball_gen_40 = [592.32,625.66,83.89,335.85,481.63,401.30,40.84,276.51]
new_ball_boot_40 = [593.32,548.71,120.78,249.98,481.63,361.49,48.96,219.02]
Det_ball_40 = [592.32,255.91,0.01,207.67,481.63,212.88,39.69,123.40]
#plot(og_ball_gen-new_ball_gen_40, label="SGM 1")
#plot!(og_ball_boot-new_ball_boot_40, label="SGM 2")

p = plot(new_ball_gen_40-Det_ball_40, label = "SGM 1",
xlabel = "Instance", ylabel = "Difference in ballast water",
title = "Ballast water difference between \n
      stochastic and deterministic solution")
plot!(p,new_ball_boot_40-Det_ball_40, label = "SGM 2")
savefig(p, plot_folder * "Ballast_water_difference_40.png")


vessel = Deterministic_problem[1].vessel
slots = Deterministic_problem[1].slots