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
    hist = histogram(weights, bins = edges, stack =true, xlabel = "Weight (t.)", ylabel = "Frequency", label = ["Type 1" "Type 2" "Type 3" "Type 4"], title = "Scenario $sc")
    return hist
end

# Plots the distrubtion of the 4 car types for the original problem
function plot_cargo_OG(problem::StowageProblem,OG::Bool = true)
    weight1 = [x.weight for x in filter(x -> x.cargo_type_id == 1, problem.cargo)]
    weight2 = [x.weight for x in filter(x -> x.cargo_type_id == 2, problem.cargo)]
    weight3 = [x.weight for x in filter(x -> x.cargo_type_id == 3, problem.cargo)]
    weight4 = [x.weight for x in filter(x -> x.cargo_type_id == 4, problem.cargo)]
    plots = []
    min_weight1 = minimum(weight1)
    max_weight1 = maximum(weight1)
    if length(weight2) > 0
        min_weight2 = minimum(weight2)
        max_weight2 = maximum(weight2)
        edges2 = min_weight2:bin_width:max_weight2
    end
    min_weight3 = minimum(weight3)
    max_weight3 = maximum(weight3)
    min_weight4 = minimum(weight4)
    max_weight4 = maximum(weight4)
    bin_width = 2 # bin width in tons
    edges1 = min_weight1:bin_width:max_weight1
    edges3 = min_weight3:bin_width:max_weight3
    edges4 = min_weight4:bin_width:max_weight4
    if length(weight2) > 0
        p = plot(
            histogram(weight1,bin = edges1, xlabel = "Weight (t.)", ylabel = "Frequency", label = "Truck"),
            histogram(weight2,bin = edges2, xlabel = "Weight (t.)", ylabel = "Frequency", label = "Car"),
            histogram(weight3, bin = edges3, xlabel = "Weight (t.)", ylabel = "Frequency", label = "Heavy Machine"),
            histogram(weight4, bin = edges4, xlabel = "Weight (t.)", ylabel = "Frequency", label = "Secu-box"),
            layout = (2,2),
            plot_title = OG ? "Cargo weight distribution for original problem" : "Cargo weight distribution for EVP"
            )
    else
        p = plot(
            histogram(weight1,bin = edges1, xlabel = "Weight (t.)", ylabel = "Frequency", label = "Truck"),
            histogram(weight3, bin = edges3, xlabel = "Weight (t.)", ylabel = "Frequency", label = "Heavy Machine"),
            histogram(weight4, bin = edges4, xlabel = "Weight (t.)", ylabel = "Frequency", label = "Secu-box"),
            layout = (2,2),
            plot_title = OG ? "Cargo weight distribution for original problem" : "Cargo weight distribution for EVP"
            )
    end
    push!(plots,p)
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
    if length(weight2)>0
        hist = histogram(weights, bins = edges, stack =true, xlabel = "Weight (t.)", ylabel = "Frequency", 
        label = ["Truck" "Car" "Heavy Machine" "Secu-box"], 
        title = OG ? "Cargo weight distribution for original problem" : "Cargo weight distribution for EVP")
    else
        hist = histogram([weight1,weight3,weight4], bins = edges, stack =true, xlabel = "Weight (t.)", 
        ylabel = "Frequency", label = ["Truck" "Heavy Machine" "Secu-box"], 
        title = OG ? "Cargo weight distribution for original problem" : "Cargo weight distribution for EVP")
    end
    push!(plots,hist)
    return plots
end

# Test normality assumption
# TODO: need guesses for means and stds for each cargo type
# NB! DOES NOT WORK CURRENTLY
# Not sure this is needed 
function test_normality(problem::StowageProblem,means, stds)
    p_values = []
    assumption = Bool[]
    weight1 = [x.weight for x in filter(x -> x.cargo_type_id == 1, problem.cargo)]
    weight2 = [x.weight for x in filter(x -> x.cargo_type_id == 2, problem.cargo)]
    weight3 = [x.weight for x in filter(x -> x.cargo_type_id == 3, problem.cargo)]
    weight4 = [x.weight for x in filter(x -> x.cargo_type_id == 4, problem.cargo)]
    for (i,weight) in enumerate([weight1,weight2,weight3,weight4])
        if length(i) > 0
            p = ExactOneSampleKSTest(weight, Normal(means(i),stds(i))) # Kolgomorov-Smirnov test
            push!(p_values,p.p_values)
            #test of the null hypothesis that the data in vector x 
            #comes from the distribution d against the alternative 
            # hypothesis that the sample is not drawn from d.
            if p < 0.05
                push!(assumption,true)
            else
                push!(assumption,false)
            end
        end
    end
end
