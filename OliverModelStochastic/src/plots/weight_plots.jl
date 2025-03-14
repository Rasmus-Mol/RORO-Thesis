# Functions to plot distribution of weight for a given scenario
function plot_cargo_weights(problem::StochasticStowageProblem,sc::Vector{Int64} = nothing)
    if isnothing(sc) # if no scenarios given plot the first 3 scenarios
        sc = [1,2,3] 
    end
    plots = []
    for i in sc
        push!(plots,plot_weights(problem,i))
    end
    return plots
end
function plot_weights(problem::StochasticStowageProblem,sc::Int64 = 1)
    weight1 = [x.weight for x in filter(x -> x.cargo_type_id == 1, problem.cargo.items[sc])]
    weight2 = [x.weight for x in filter(x -> x.cargo_type_id == 2, problem.cargo.items[sc])]
    weight3 = [x.weight for x in filter(x -> x.cargo_type_id == 3, problem.cargo.items[sc])]
    weight4 = [x.weight for x in filter(x -> x.cargo_type_id == 4, problem.cargo.items[sc])]
    weights = [weight1,weight2,weight3,weight4]
    min_weight = []
    max_weight = []
    for i in 1:length(weights)
        if length(weights[i]) > 0
            push!(min_weight,minimum(weights[i]))
            push!(max_weight,maximum(weights[i]))
        end
    end
    min_weight = minimum(min_weight)
    max_weight = maximum(max_weight)
    bin_width = 2 # bin width in tons
    edges = min_weight:bin_width:max_weight
    hist = histogram(weights, bins = edges, stack =true, xlabel = "Weight", ylabel = "Frequency", label = ["Type 1" "Type 2" "Type 3" "Type 4"], title = "Cargo weight distribution for scenario $sc")
    return hist
end

# Plots the distrubtion of the 4 car types for the original problem
function plot_cargo_OG(problem::StowageProblem)
    weight1 = [x.weight for x in filter(x -> x.cargo_type_id == 1, problem.cargo)]
    #weight2 = [x.weight for x in filter(x -> x.cargo_type_id == 2, problem.cargo)]
    weight3 = [x.weight for x in filter(x -> x.cargo_type_id == 3, problem.cargo)]
    weight4 = [x.weight for x in filter(x -> x.cargo_type_id == 4, problem.cargo)]
    #weights = [weight1,weight2,weight3,weight4]
    min_weight1 = minimum(weight1)
    max_weight1 = maximum(weight1)
    #min_weight2 = minimum(weight2)
    #max_weight2 = maximum(weight2)
    min_weight3 = minimum(weight3)
    max_weight3 = maximum(weight3)
    min_weight4 = minimum(weight4)
    max_weight4 = maximum(weight4)
    bin_width = 2 # bin width in tons
    edges1 = min_weight1:bin_width:max_weight1
    #edges2 = min_weight2:bin_width:max_weight2
    edges3 = min_weight3:bin_width:max_weight3
    edges4 = min_weight4:bin_width:max_weight4
    p = plot(
        histogram(weight1,bin = edges1, xlabel = "Weight", ylabel = "Frequency", label = "Type 1"),
        #histogram(weight2,bin = edges2, xlabel = "Weight", ylabel = "Frequency", label = "Type 2"),
        histogram(weight3, bin = edges3, xlabel = "Weight", ylabel = "Frequency", label = "Type 3"),
        histogram(weight4, bin = edges4, xlabel = "Weight", ylabel = "Frequency", label = "Type 4"),
        layout = (2,2),
        plot_title = "Cargo weight distribution for original problem"
    )
    return p
end

