# Estimate effect of scenario reduction
include("packages_and_files.jl")

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
# add biased noise
# Create problem using few scenarios and no scenario reduction
# Create problem using many scenarios and using scenario reduction
# Estimate which is closest to Deterministic problem ?
for i in 1:8 # number of instances
    problem_det_noise = add_biased_noise_to_test_instance(Deterministic_problem[i])
    sto_no_reduction[i] = create_stochastic_problem(problem_det_noise, sc, length(Deterministic_problem[i].cargo), [], Bootstrap_bookedweight_quantile)  
    for j in 1:length(scenario_reduction)
        pro = create_stochastic_problem(problem_det_noise, scenario_reduction[j], length(Deterministic_problem[i].cargo), [], Bootstrap_bookedweight_quantile)
        sto_reduction[i,j] = scenario_reduced(pro, sc, scenario_reduction_clustering, 60)
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
p = plot(avg_dist, label = "No reduction", title = "Average distance from original problem", xlabel = "Instance", ylabel = "Distance")
for j in 1:length(scenario_reduction)
    plot!(p, avg_dist_reduction[:,j], label = "Initialscenarios $(scenario_reduction[j])")
end
display(p)
# distance to original problem
n = length.(dist_no_reduction)   # number of samples per instance
means = mean.(dist_no_reduction)
stds = std.(dist_no_reduction)
stderr = stds ./ sqrt.(n)
# 95% confidence interval z-value
z = 1.96
lower = means .- z .* stderr
upper = means .+ z .* stderr
yerr = z .* stderr
# Plot with error bars
plot(
    means,
    yerror = yerr,
    #ribbon = (means .- lower, upper .- means),
    xlabel = "Instance",
    ylabel = "Mean distance from original problem",
    title = "No scenario reduction",
    label = ""
    #label = "Mean ± CI",
    #legend = :topright
)
for i in 1:length(scenario_reduction)
    n = length.(dist_reduction[:,i])   # number of samples per instance
    means = mean.(dist_reduction[:,i])
    stds = std.(dist_reduction[:,i])
    stderr = stds ./ sqrt.(n)
    # 95% confidence interval z-value
    z = 1.96
    lower = means .- z .* stderr
    upper = means .+ z .* stderr
    yerr = z .* stderr
    # Plot with error bars
    p = plot(
        means,
        yerror = yerr,
        #ribbon = (means .- lower, upper .- means),
        xlabel = "Instance",
        ylabel = "Mean distance from original problem",
        title = "Scenario reduction - initial scenarios $(scenario_reduction[i])",
        label = ""
        #label = "Mean ± CI",
        #legend = :topright
    )
    display(p)
end

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
