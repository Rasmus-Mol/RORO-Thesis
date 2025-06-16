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
sol_evp = Array{Any}(nothing, 5,6,7,8)
nested = [[[ [sol_evp[1,j,k,l] for l in 1:size(sol_evp,4)]
for k in 1:size(sol_evp,3)]
for j in 1:size(sol_evp,2)] for i in 1:size(sol_evp,1)]

open("test_EVP.json", "w") do file
    JSON.print(file, nested, 4)
end
function get_obj_of_evp(filename::String)
    open(filename*".json", "r") do file
        obj = JSON3.read(read(file, String))#, Vector{Vector{Vector{Vector{Float64}}}})
        return obj
    end
end
as = get_obj_of_evp("test_EVP")
ass = to_4d_array(as)

length(as)
function to_4d_array(nested)
    dims = (
        length(nested),
        length(nested[1]),
        length(nested[1][1]),
        length(nested[1][1][1])
    )
    A = Array{Any, 4}(undef, dims...)
    for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3], l in 1:dims[4]
        A[i, j, k, l] = nested[i][j][k][l]
    end
    return A
end
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