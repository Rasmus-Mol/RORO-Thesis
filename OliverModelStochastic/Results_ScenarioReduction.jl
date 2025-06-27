# Estimate effect of scenario reduction
include("packages_and_files.jl")
plot_folder = "Plots/Results/BiasedNoise"
# Create folder for plots
if !isdir(plot_folder)
    mkpath(plot_folder)
end

Deterministic_problem = Array{Any}(undef, 8)
problemname1, problemname3 = "finlandia", "hazardous"
for i in 1:8
    Deterministic_problem[i] = load_data(problemname1, Finlandia_test[i], problemname3)
end
################################
vessel = Deterministic_problem[1].vessel
for i in 1:length(vessel.ballast_tanks)
    println("Max vcg: ", vessel.ballast_tanks[i].max_vcg)
    println("My slope: ", round((vessel.ballast_tanks[i].max_vcg - vessel.ballast_tanks[i].min_vcg) / vessel.ballast_tanks[i].max_vol,digits = 2)
    )
end
vessel.ballast_tanks[1].lcg
vessel.ballast_tanks[1].tcg
vessel.ballast_tanks[2].min_vcg
vessel.ballast_tanks[2].max_vcg
vessel.ballast_tanks[6].max_vol
length(vessel.ballast_tanks)
(vessel.ballast_tanks[2].max_vcg-vessel.ballast_tanks[2].min_vcg)/vessel.ballast_tanks[6].max_vol
HPC_folder = "Finlandia_mixed_light_60_28_05_13/"
Deterministic_Solution = get_solution_deterministic("Finlandia_deterministic",
"Deterministic_Solution",HPC_folder)
#####################
# What issues does an empty ship have.
test_problem_name = Finlandia_test[1]
HPC_folder_load = HPC_folders[1]
HPC_folder_save = HPC_folder_load
# Load solution
det_sol = get_solution_deterministic("Finlandia_deterministic",
        "Deterministic_Solution", HPC_folder_load)
det_pro = load_data(problemname1, test_problem_name, problemname3)
cs_empty = zeros(length(det_pro.cargo),length(det_pro.slots))
model_empty_no_slack = second_stage_model(cs_empty, det_pro)
set_time_limit_sec(model_empty_no_slack, 60 * 60) # 60 minutes
set_silent(model_empty_no_slack)
optimize!(model_empty_no_slack)
termination_status(model_empty_no_slack)
sum(value.(model_empty_no_slack[:ballast_volume]))
sum(getfield.(det_pro.vessel.ballast_tanks, :max_vol))
# Different slope
model_empty_no_slack_old_slope = second_stage_model_different_slope(cs_empty, det_pro)
set_time_limit_sec(model_empty_no_slack_old_slope, 60 * 60) # 60 minutes
set_silent(model_empty_no_slack_old_slope)
optimize!(model_empty_no_slack_old_slope)
termination_status(model_empty_no_slack_old_slope)
sum(value.(model_empty_no_slack_old_slope[:ballast_volume]))
# same ballast water
println(value.(model_empty_no_slack_old_slope[:ballast_volume]))
println(value.(model_empty_no_slack[:ballast_volume]))
# vcg_ballast
my_slop = calculate_vcg_slopes(det_pro.vessel)
println(value.(model_empty_no_slack_old_slope[:vcg_ballast]))
bal = collect(value.(model_empty_no_slack[:ballast_volume]))
println(sum(my_slop .* bal))
#################################

# Scenario reduction
scenarios = [10]
sc = scenarios[1]
scenario_reduction = [100,300,500]
sto_no_reduction = Array{Any}(undef, 8, length(scenario_reduction))
sto_reduction = Array{Any}(undef, 8, length(scenario_reduction))
pro_reduction = Array{Any}(undef, 8, length(scenario_reduction))
# add biased noise
# Create problem using few scenarios and no scenario reduction
# Create problem using many scenarios and using scenario reduction
# Estimate which is closest to Deterministic problem ?
for i in 1:8 # number of instances
    println("it: ", i)
    problem_det_noise = add_biased_noise_to_test_instance(Deterministic_problem[i])
    sto_no_reduction[i] = create_stochastic_problem(problem_det_noise, sc, length(Deterministic_problem[i].cargo), [], Bootstrap_bookedweight_quantile)  
    for j in 1:length(scenario_reduction)
        pro_reduction[i,j] = create_stochastic_problem(problem_det_noise, scenario_reduction[j], length(Deterministic_problem[i].cargo), [], Bootstrap_bookedweight_quantile)
        sto_reduction[i,j] = scenario_reduced(pro_reduction[i,j] , sc, scenario_reduction_clustering, 60)
    end
end

dist_no_reduction = Array{Any}(undef, 8)
dist_reduction = Array{Any}(undef, 8, length(scenario_reduction))
for i in 1:8 # number of instances
    dist_no_reduction[i] = distance_from_original_problem(sto_no_reduction[i], Deterministic_problem[i])
    for j in 1:length(scenario_reduction)
        dist_reduction[i,j] = distance_from_original_problem(sto_reduction[i,j], Deterministic_problem[i])
    end
end
println(dist_no_reduction[1])
println(dist_reduction[1,1])

for i in 1:length(scenario_reduction)
    row_names = ["Instance $(j)" for j in 1:8]
    header = ["Scenario $(j)" for j in 1:length(dist_reduction[1,i])]
    temp =reshape(vcat(dist_reduction[:,i]...), length(dist_reduction[:,i][1]), length(dist_reduction[:,i]))'
    temp1 = DataFrame(temp, :auto)
    rename!(temp1, Symbol.(header))
    insertcols!(temp1,1,:Instance =>row_names)
    pretty_table(temp1)
end

avg_dist = [mean(dist_no_reduction[i]) for i in 1:8]
avg_dist_reduction = [mean(dist_reduction[i,j]) for i in 1:8, j in 1:length(scenario_reduction)]
p = plot(avg_dist, label = "No reduction", title = "Average distance to deterministic instance weight", xlabel = "Instance", ylabel = "Distance")
for j in 1:length(scenario_reduction)
    plot!(p, avg_dist_reduction[:,j], label = "Initialscenarios $(scenario_reduction[j])")
end
display(p)
savefig(p, plot_folder*"/Average_distance_to_deterministic_instance_weight.png")
# distance to original problem
n = length.(dist_no_reduction)   # number of samples per instance
means = mean.(dist_no_reduction)
stds = std.(dist_no_reduction)
stderrs = stds ./ sqrt.(n)
# 95% confidence interval z-value
z = 1.96
lower = means .- z .* stderrs
upper = means .+ z .* stderrs
yerr = z .* stderrs
# Plot with error bars
p1 = plot(
    means,
    yerror = yerr,
    #ribbon = (means .- lower, upper .- means),
    xlabel = "Instance",
    ylabel = "Mean distance",
    title = "No scenario reduction",
    label = ""
    #label = "Mean ± CI",
    #legend = :topright
)
savefig(p1, plot_folder*"/Mean_distance_No_scenario_reduction.png")
for i in 1:length(scenario_reduction)
    n = length.(dist_reduction[:,i])   # number of samples per instance
    means = mean.(dist_reduction[:,i])
    stds = std.(dist_reduction[:,i])
    stderrs = stds ./ sqrt.(n)
    # 95% confidence interval z-value
    z = 1.96
    lower = means .- z .* stderrs
    upper = means .+ z .* stderrs
    yerr = z .* stderrs
    # Plot with error bars
    p = plot(
        means,
        yerror = yerr,
        #ribbon = (means .- lower, upper .- means),
        xlabel = "Instance",
        ylabel = "Mean distance",
        title = "Scenario reduction - initial scenarios $(scenario_reduction[i])",
        label = ""
        #label = "Mean ± CI",
        #legend = :topright
    )
    display(p)
    savefig(p, plot_folder*"/Mean_distance_scenario_reduction_sc_$(scenario_reduction[i]).png")
end
# Find silhouette scores
# for 1 instance
temp = silhouette_score(pro_reduction[1,1].cargo ,10)
clusters = [5,10,15]
s_means = zeros(length(clusters),8, length(scenario_reduction))
for i in 1:8
    println("Instance: ", i)
    for j in 1:length(scenario_reduction)
        for l in 1:length(clusters) # different number of clusters
            s_means[l,i,j] = mean(silhouette_score(pro_reduction[i,j].cargo ,clusters[l]))
        end
    end
end
for i in 1:8
    println("##############")
    println("instance: ", i)
    println(round.(s_means[2,i,:], digits=4))
end
s_means_instances = round.(mean(s_means, dims = 2), digits = 4)
println(s_means_instances[:,1,:])
s_means[end,:,:]
CargoC = pro_reduction[5,1].cargo
sc = length(CargoC.items)
cost = zeros(sc, sc)
for i in 1:sc
    for j in 1:sc
        cost[i,j] = sum(abs.(getfield.(CargoC.items[i],:weight) .- getfield.(CargoC.items[j],:weight)))
    end
end
# Perform hierarchical clustering
result = hclust(cost, linkage = :complete) #:average)  # or :single, :complete, :ward, etc.
# Cut the dendrogram into sc clusters
labels = cutree(result, k = 10)
scores = silhouette_score(pro_reduction[5,1].cargo ,10)
# plot silhouette scores
plot_silhouette(scores, labels)
plot_silhouette2(scores, labels)
mean(scores)
scores


###############
test_problem_name = Finlandia_test[1]
problem_det = load_data("finlandia", test_problem_name, "hazardous")
problem_det_noise = add_biased_noise_to_test_instance(problem_det)
pro = create_stochastic_problem(problem_det_noise, 100, length(problem_det.cargo), [], Bootstrap_bookedweight_quantile)
pro = scenario_reduced(pro, 10, scenario_reduction_clustering, 60)
mo = create_model_stochastic(pro)
    set_silent(mo) # removes terminal output
    set_time_limit_sec(mo, 5*60)
    optimize!(mo)
pro = create_stochastic_problem_cars_known(problem_det_noise, 10)
pro.cargo.items[1].items[50].weight
pro.cargo.items[2].items[50].weight
for i in 1:10
    println(pro.cargo.items[i].items[50].weight)
end
problem_det.cargo.items[50].weight
###############

function distance_from_original_problem(sto_pro::StochasticStowageProblem, det_pro::StowageProblem)
    # Calculate distance between the cargoes in the problem and the original problem
    sc = sto_pro.scenarios
    dist = zeros(sc)
    det_truckssecu = filter(x-> x.cargo_type_id == 1 || x.cargo_type_id == 4, det_pro.cargo.items)
    for i in 1:sc
        # filter such that it is only trucks and Secu
        sto_truckssecu = filter(x-> x.cargo_type_id == 1 || x.cargo_type_id == 4, sto_pro.cargo.items[i].items)
        dist[i] = round(sum(abs.(getfield.(sto_truckssecu,:weight) .- getfield.(det_truckssecu,:weight))),digits = 2)
    end
    return dist
end

function silhouette_score(CargoC::CargoCollectionScenarios, no_clusters)
    sc = length(CargoC.items)
    cost = zeros(sc, sc)
    for i in 1:sc
        for j in 1:sc
            cost[i,j] = sum(abs.(getfield.(CargoC.items[i],:weight) .- getfield.(CargoC.items[j],:weight)))
        end
    end
    # Perform hierarchical clustering
    result = hclust(cost, linkage = :complete) #:average)  # or :single, :complete, :ward, etc.
    # Cut the dendrogram into sc clusters
    labels = cutree(result, k = no_clusters)
    # Find silhouette scores
    n = size(cost, 1)
    s = zeros(n)

    for i in 1:n
        same_cluster = findall(j -> labels[j] == labels[i] && j != i, 1:n)
        other_clusters = setdiff(unique(labels), [labels[i]])

        # a(i): avg distance to same cluster
        a_i = isempty(same_cluster) ? 0.0 : mean(cost[i, same_cluster])

        # b(i): min avg distance to other clusters
        b_i = minimum([
            mean(cost[i, findall(j -> labels[j] == c, 1:n)])
            for c in other_clusters
        ])

        s[i] = if a_i == b_i
            0.0
        elseif a_i < b_i
            1 - a_i / b_i
        else
            b_i / a_i - 1
        end
        if a_i == 0.0
            s[i] = 0
        end
    end
    return s
end

function plot_silhouette(scores::Vector{Float64}, labels::Vector{Int})
    n = length(scores)
    cluster_ids = sort(unique(labels))
    cluster_scores = [scores[findall(labels .== c)] for c in cluster_ids]
    
    # Sort each cluster's scores (optional, for prettier plot)
    sorted_scores = [sort(s, rev=true) for s in cluster_scores]

    # Prepare bar data
    y = Float64[]
    cluster_colors = Int[]
    x = Int[]
    tick_positions = Int[]
    tick_labels = String[]
    pos = 1

    for (i, scs) in enumerate(sorted_scores)
        append!(y, scs)
        append!(cluster_colors, fill(i, length(scs)))
        append!(x, pos:pos+length(scs)-1)
        push!(tick_positions, pos + length(scs) ÷ 2)
        push!(tick_labels, "Cluster $(cluster_ids[i])")
        pos += length(scs)
    end

    bar(x, y, group=cluster_colors, legend=false, xlabel="Data points", ylabel="Silhouette Score",
        title="Silhouette Scores by Cluster", xticks=(tick_positions, tick_labels))
end

function plot_silhouette2(scores::Vector{Float64}, labels::Vector{Int})
    n = length(scores)
    cluster_ids = sort(unique(labels))
    cluster_scores = [scores[findall(labels .== c)] for c in cluster_ids]
    
    # Sort scores within each cluster (optional for aesthetics)
    sorted_scores = [sort(s, rev=true) for s in cluster_scores]

    y = Float64[]
    cluster_colors = Int[]
    x = Int[]
    tick_positions = Int[]
    tick_labels = String[]
    pos = 1

    for (i, scs) in enumerate(sorted_scores)
        append!(y, scs)
        append!(cluster_colors, fill(i, length(scs)))
        append!(x, pos:pos+length(scs)-1)
        push!(tick_positions, pos + length(scs) ÷ 2)
        push!(tick_labels, "Cluster $(cluster_ids[i])")
        pos += length(scs)
    end

    bar(x, y, group=cluster_colors, legend=false,
        xlabel="Data Points", ylabel="Silhouette Score",
        title="Silhouette Scores by Cluster",
        xticks=(tick_positions, tick_labels),
        xrotation=45)  # << Rotate x-axis labels by 45 degrees
end

########################################################################
# Load results
########################################################################
HPC_folders_det = [
    "Finlandia_mixed_light_60_28_05_13",
    "Finlandia_mixed_light_100_28_05_13",
    "Finlandia_mixed_heavy_60_28_05_14",
    "Finlandia_mixed_heavy_100_28_05_15",
    "Finlandia_no_cars_light_60_28_05_16",
    "Finlandia_no_cars_light_100_28_05_16",
    "Finlandia_no_cars_heavy_60_28_05_19",
    "Finlandia_no_cars_heavy_100_28_05_19",
]
# load deterministic results
problemname1, problemname3 = "finlandia", "hazardous"
Deterministic_Solution = Array{Any}(undef, length(HPC_folders_det))
Deterministic_problem = Array{Any}(undef, length(HPC_folders_det))
for i in 1:length(HPC_folders_det)
    HPC_folder = HPC_folders_det[i]
    Deterministic_Solution[i] = get_solution_deterministic("Finlandia_deterministic",
    "Deterministic_Solution", HPC_folder)
    test_instance = Finlandia_test[i]
    Deterministic_problem[i] = load_data(problemname1, test_instance, problemname3)
end
# load stochastic solutions
function Get_EVP_results_2D(HPC_folder::String, filename::String)
    temp = open(joinpath("Results",HPC_folder,filename*".json"), "r") do file
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
HPC_folders = [
    "Finlandia_mixed_light_60_sc_100_biasednoise_16_06_22",
    "Finlandia_mixed_light_100_sc_100_biasednoise_16_06_22",
    "Finlandia_mixed_heavy_60_sc_100_biasednoise_16_06_22",
    "Finlandia_mixed_heavy_100_sc_100_biasednoise_16_06_22",
    "Finlandia_no_cars_light_60_sc_100_biasednoise_16_06_22",
    "Finlandia_no_cars_light_100_sc_100_biasednoise_17_06_09",
    "Finlandia_no_cars_heavy_60_sc_100_biasednoise_16_06_22",
    "Finlandia_no_cars_heavy_100_sc_100_biasednoise_16_06_22",
]

scenarios = [10]
repetitions = 5 
EVP_boot = Array{Any}(undef, length(HPC_folders),repetitions,)
EVP_boot_matrix = Array{Any}(undef, length(HPC_folders), 4)
Stochastic_boot = Array{Any}(undef, length(HPC_folders), repetitions)
Stochastic_boot_fitted = Array{Any}(undef, length(HPC_folders), repetitions)

boot_fitted_slacked = Array{Any}(nothing, length(HPC_folders), repetitions)
slack_boot = Array{Any}(nothing, length(HPC_folders), repetitions)

for i in 1:length(HPC_folders)
    test_instance = Finlandia_test[i]
    HPC_folder = HPC_folders[i]
    EVP_boot_matrix[i,1] = Get_EVP_results_2D(HPC_folders[i], "Objective_values")
    EVP_boot_matrix[i,2] = Get_EVP_results_2D(HPC_folders[i], "Ballast_used")
    EVP_boot_matrix[i,3] = Get_EVP_results_2D(HPC_folders[i], "Cargo_loaded")
    EVP_boot_matrix[i,4] = Get_EVP_results_2D(HPC_folders[i], "Inf_index")
    for j in 1:repetitions
        foldername = "EVP_Bootstrap1_rep$(j)_sc10_unknown$(Deterministic_Solution[i].n_cargo_total)_time3600"
        filename = "EVP_Solution"
        EVP_boot[i,j] = get_solution_deterministic(foldername,filename,HPC_folder)
        foldername = "Stochastic_Bootstrap1_rep$(j)_sc$(10)_unknown$(Deterministic_Solution[i].n_cargo_total)_time3600"
        filename = "Stochastic_Solution"
        Stochastic_boot[i,j] = get_solution_stochastic(foldername,
        filename, HPC_folder)
        Stochastic_boot_fitted[i,j] = get_solution_deterministic(foldername,
                    "Fitted_Solution", HPC_folder)
        # Load slacked
        if Stochastic_boot_fitted[i,j].status != "OPTIMAL"
            boot_fitted_slacked[i,j] = get_solution_deterministic(foldername,
            "Fitted_Solution_slacked", HPC_folder)
            slack_boot[i,j] = get_slack(foldername, "Fitted_Solution",HPC_folder)
        end
    end
end
# EVP_matrix: obj,ballast, cargo, infeasibility
# no infeasibility - all objective values are correct
corrected_obj = zeros(length(HPC_folders), repetitions, 10)
base_line_obj = [-500000,-900000,-500000,-900000,-200000,-300000,-200000,-300000 ]
EVP_inf_count = zeros(length(HPC_folders),repetitions)
# Some test were maybe infeasible, finding objective value for all tests
# If test were feasible, objective value is the same
for i in 1:length(HPC_folders)
    for j in 1:repetitions
        for l in 1:10
            corrected_obj[i,j,l] = get_obj_val(Float64(EVP_boot_matrix[i,3][j,l]), 
            Deterministic_problem[i], EVP_boot_matrix[i,2][j,l])
            if abs(corrected_obj[i,j,l] - EVP_boot_matrix[i,1][j,l]) > 0.01
                EVP_inf_count[i,j] += 1
            end
            ## test was infeasible
            #if EVP_boot_matrix[i,1][j,l] > base_line_obj[i]
            #    corrected_obj
            #else
            #end
        end
    end
end
EVP_boot_matrix[1,3]
isnothing(idx_inf[1,1])
println(EVP_boot_matrix[8,1])
EVP_boot_matrix[8,1]
corrected_obj[8,:,:]
EVP_inf_count
# EEV
EEV = dropdims(mean(corrected_obj, dims = 3); dims = 3)
#EEV_ballast_water = [dropdims(mean(EVP_boot_matrix[i,2], dims = 2);dims=2) for i in 1:length(HPC_folders)]
sto_obj = zeros(length(HPC_folders), repetitions)
sto_fitted_obj = zeros(length(HPC_folders), repetitions)
#sto_ballast_water
for i in 1:length(HPC_folders)
    for j in 1:repetitions
        sto_obj[i,j] = Stochastic_boot[i,j].objective
        if Stochastic_boot_fitted[i,j].objective == Inf # infeasible
            println("INF: ", i, j)
            sto_fitted_obj[i,j] = get_obj_val(
                Float64(boot_fitted_slacked[i,j].n_cargo_loaded),
                Deterministic_problem[i],
                boot_fitted_slacked[i,j].ballast_weight
            )
        else # feasible
            sto_fitted_obj[i,j] = Stochastic_boot_fitted[i,j].objective
        end
    end
end
(Stochastic_boot_fitted[4,2].ballast_weight+Stochastic_boot_fitted[4,3].ballast_weight+
Stochastic_boot_fitted[4,4].ballast_weight+Stochastic_boot_fitted[4,5].ballast_weight)/4
Stochastic_boot_fitted[8,5].objective 
# Checking unique solutions
for i in 1:length(HPC_folders)
    println(length(unique(getfield.(Stochastic_boot[i,:],:cs))))
end

# comparing ballast water used to deterministic
ballast_w = Array{Any}(undef, length(HPC_folders))
for i in 1:length(HPC_folders)
    temp = getfield.(Stochastic_boot_fitted[i,:], :ballast_weight)
    temp = filter(x -> x != 0.0, temp) # remove zero ballast water
    if length(temp)>0
        ballast_w[i] = mean(temp) - Deterministic_Solution[i].ballast_weight
    else
        ballast_w[i] = [0.0]
    end
    #for j in 1:repetitions
    #    if Stochastic_boot_fitted[i,j].ballast_weight == 0.0
    #        deleteat!(ballast_w[i], j)
    #        j -= 1
    #    end
    #end
end
x = [i for i in 1:length(HPC_folders) if ballast_w[i] != [0.0]]
p = plot(x,round.(filter(x-> x!= [0.0], ballast_w),digits = 3), label = "", 
    title = "Mean ballast water difference between \n
    stochastic and deterministic solution",
    xlabel = "Instance", ylabel = "Difference in ballast water")
savefig(p, plot_folder*"/Ballast_water_difference.png")
# VVS
VSS = EEV .- sto_obj
#VSS_2 = EEV .- sto_fitted_obj
getfield.(Stochastic_boot_fitted[7,:],:ballast_weight)
getfield.(Stochastic_boot[2,:],:gap)
getfield.(Stochastic_boot[4,:],:gap)
Stochastic_boot[8,1].n_cargo_loaded
Stochastic_boot[8,3].n_cargo_loaded
VSS_ballast_water = 
EVP_boot_matrix[8,1]
EEV[8,1]
VSS
sto_obj
# EVPI
# If value is above zero -> it is because sto_sol was infeasible
# and infeasible placement was used and constraint was slacked and 
# solution found was better
# Dropping the infeasible solutions
# (4,1) & (8,:)
EVPI = sto_fitted_obj .- getfield.(Deterministic_Solution[:], :objective)
EVPI
EVPI_2 = sto_obj .- getfield.(Deterministic_Solution[:], :objective)
EVPI_2
# Plot EVPI
p1 = boxplot([EVPI[1,:], EVPI[3,:], EVPI[5,:], EVPI[6,:],EVPI[7,:]]; 
        #labels = ["Instance 1", "Instance 3", "Instance 5 C",
        #"Instance 6", "Instance 7"],    # custom labels for each box
        labels = "",
        xticks = (1:5, ["1", "3", "5", "6", "7"]),
        #color = [:teal, :orange, :purple],            # fill colors for the boxes
        xlabel = "Instance", 
        ylabel = "EVPI",           # axis labels
        title = "Expected Value of Perfect Information")
savefig(p1,plot_folder*"/EVPI_1_3_5_6_7.png")
p2 = boxplot([EVPI[8,:]]; 
        #labels = ["Instance 1", "Instance 3", "Instance 5 C",
        #"Instance 6", "Instance 7"],    # custom labels for each box
        labels = "",
        xticks = (1:1, ["8"]),
        #color = [:teal, :orange, :purple],            # fill colors for the boxes
        xlabel = "Instance", 
        ylabel = "EVPI",           # axis labels
        title = "Expected Value of Perfect Information")
savefig(p2,plot_folder*"/EVPI_8.png")
p3 = boxplot([EVPI[8,2:5]]; 
        #labels = ["Instance 1", "Instance 3", "Instance 5 C",
        #"Instance 6", "Instance 7"],    # custom labels for each box
        labels = "",
        xticks = (1:1, ["8"]),
        #color = [:teal, :orange, :purple],            # fill colors for the boxes
        xlabel = "Instance", 
        ylabel = "EVPI",           # axis labels
        title = "Expected Value of Perfect Information")
savefig(p3,plot_folder*"/EVPI_8_partly.png")
p4 = boxplot([EVPI[2,:], EVPI[4,2:5]]; 
        #labels = ["Instance 1", "Instance 3", "Instance 5 C",
        #"Instance 6", "Instance 7"],    # custom labels for each box
        labels = "",
        xticks = (1:2, ["2", "4"]),
        #color = [:teal, :orange, :purple],            # fill colors for the boxes
        xlabel = "Instance", 
        ylabel = "EVPI",           # axis labels
        title = "Expected Value of Perfect Information")
savefig(p4,plot_folder*"/EVPI_2_4.png")

######################################################
# plots VVS

# Create a basic boxplot for the three groups
p1 = boxplot([VSS[1,:], VSS[3,:], VSS[5,:], VSS[6,:],VSS[7,:]]; 
        #labels = ["Instance 1", "Instance 3", "Instance 5 C",
        #"Instance 6", "Instance 7"],    # custom labels for each box
        labels = "",
        xticks = (1:5, ["1", "3", "5", "6", "7"]),
        #color = [:teal, :orange, :purple],            # fill colors for the boxes
        xlabel = "Instance", 
        ylabel = "VSS",           # axis labels
        title = "Value of Stochastic Solution")
savefig(p1,plot_folder*"/VSS_1_3_5_6_7.png")
p2 = boxplot([VSS[8,:]]; 
        #labels = ["Instance 1", "Instance 3", "Instance 5 C",
        #"Instance 6", "Instance 7"],    # custom labels for each box
        labels = "",
        xticks = (1, ["8"]),
        #color = [:teal, :orange, :purple],            # fill colors for the boxes
        xlabel = "Instance", 
        ylabel = "VSS",           # axis labels
        title = "Value of Stochastic Solution")
savefig(p2,plot_folder*"/VSS_8.png")
p3 = boxplot([VSS[8,2:5]]; 
        #labels = ["Instance 1", "Instance 3", "Instance 5 C",
        #"Instance 6", "Instance 7"],    # custom labels for each box
        labels = "",
        xticks = (1, ["8"]),
        #color = [:teal, :orange, :purple],            # fill colors for the boxes
        xlabel = "Instance", 
        ylabel = "VSS",           # axis labels
        title = "Value of Stochastic Solution")
savefig(p3,plot_folder*"/VSS_8_partly.png")
p4 = boxplot([VSS[2,:], VSS[4,:]]; 
        #labels = ["Instance 1", "Instance 3", "Instance 5 C",
        #"Instance 6", "Instance 7"],    # custom labels for each box
        labels = "",
        xticks = (1:2, ["2", "4"]),
        #color = [:teal, :orange, :purple],            # fill colors for the boxes
        xlabel = "Instance", 
        ylabel = "VSS",           # axis labels
        title = "Value of Stochastic Solution")
savefig(p4,plot_folder*"/VSS_2_4.png")

# results for table
for i in 1:length(HPC_folders)
    temp = zeros(repetitions)
    feasibilit_count = 0
    inf_reason = zeros(repetitions)
    for j in 1:repetitions
        #temp[j] = Stochastic_boot[i,j].gap
        #temp[j] = Stochastic_boot[i,j].n_cargo_loaded
        if isnothing(slack_boot[i,j])
            feasibilit_count +=1
            temp[j] = Stochastic_boot_fitted[i,j].ballast_weight
        else
            temp[j] = boot_fitted_slacked[i,j].ballast_weight
            if sum(slack_boot[i,j][1])>0
                inf_reason[j] = 1
            else
                inf_reason[j] = 2
            end
        end
    end
    println("################")
    #println("Instance $(i) - mean gap: ", round(mean(temp),digits = 4))
    #println("Instance $(i) - gap std: ", round(std(temp),digits = 2))
    #println("Instance $(i) - mean cargo loaded: ", round(mean(temp),digits = 4))
    #if feasibilit_count > 0
        temp2 = filter(x -> x != 0, temp)
        println("mean ballast water: ", round(mean(temp2),digits = 2))
    #end
    println("Instance $(i) - inf reasons: ", inf_reason)
end



# Finds objective value
base_line_obj = [-500000,-900000,-500000,-900000,-200000,-300000,-200000,-300000 ]
function get_obj_val(cargo_loaded::Float64, det_pro::StowageProblem, ballast_used::Float64)
    vessel = det_pro.vessel
    n_ballast_tanks = length(vessel.ballast_tanks)
    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
	#haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
	#cost = [CSC for c in cargo]
	#cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
	#area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]
    return ballast_used - CSC*cargo_loaded
end
