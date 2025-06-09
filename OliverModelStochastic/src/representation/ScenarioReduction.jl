# Scenario reduction
function scenario_reduced(problem::StochasticStowageProblem, sc::Int64, 
    reduction_method::Function = scenario_reduction_naive, timelimit::Int64 = 60)
    cargoC = problem.cargo
    pr = problem.probability
    cargoc, probability = reduction_method(cargoC, pr, sc, timelimit)
    # create new problem
    #println("prob:", probability)
    new_problem = StochasticStowageProblem(
        vessel = problem.vessel,
        slots = problem.slots,
        cargo = cargoc,
        unknown_weights = problem.unknown_weights,
        known_weights = problem.known_weights,
        scenarios = sc,
        probability = probability,
        name = problem.name,
        timestamp = problem.timestamp
    )
    return new_problem
end
# Contract until the number of scenarios is equal to sc
function scenario_reduction_naive(cargoC::CargoCollectionScenarios,
    probability::Vector{Float64}, sc::Int64, time_limit::Int64 = 60)
    sc_old = length(cargoC.items)
    if sc_old <= sc
        throw("No reduction possible")
        return cargoC, probability
    end
    scenarios = sc_old
    while scenarios != sc
        x , obj_val = Find_one_reduction_naive(cargoC, scenarios)
        if x == 0
            throw("No reduction possible")
            return cargoC
        end
        # Contract scenarios
        #println("Reducing")
        #println("")
        idx = findall(t -> t==1, x)
        sc_idx = [idx[1][1], idx[1][2]]
        # shuffle the order we contract, i.e. either i into j or j into i
        # Dont know if this is a good idea
        shuffle!(sc_idx)
        #println(scenarios)
        cargoC, probability = contract_scenarios_naive(cargoC, probability, sc_idx, scenarios)
        #println("Number of scenarios: ", length(cargoC.items))
        #println("Probabilities: ", probability)
        #println("Total weights: ", cargoC.items.total_weight)
        scenarios -= 1
        #println("stop criteria: ", scenarios != sc)
    end
    return cargoC, probability
end

function contract_scenarios_naive(cargoC::CargoCollectionScenarios, pr::Vector{Float64},
    sc, sc_old::Int64)
    # Create a new cargo collection with the reduced number of scenarios
    cargo_scenarios = Vector{CargoCollection}()
    cargoes = Vector{Cargo}()
    probability = zeros(sc_old) # array with probabilities
    for i in 1:sc_old
        #println("i:" ,i)
        #println("pr: ", pr[i])
        #println("probability: ", probability[i])
        if i == sc[1] # contract
            for j in 1:length(cargoC.items[i].items)
                # create new cargo
                current_cargo = cargoC.items[i].items[j]
                new_cargo = Cargo(
                    id=current_cargo.id,
                    cargo_type_id=current_cargo.cargo_type_id,  # Use the id from type_info
                    weight= (current_cargo.weight + cargoC.items[sc[2]].items[j].weight)/2,
                    loading_port=current_cargo.loading_port,
                    discharge_port=current_cargo.discharge_port,
                    priority=current_cargo.priority,
                    requires_lashing=current_cargo.requires_lashing,
                    requires_ventilation=current_cargo.requires_ventilation,
                    hazardous=current_cargo.hazardous,
                    refers=current_cargo.refers
                )
                push!(cargoes, new_cargo)
            end
            probability[i] = pr[sc[1]] + pr[sc[2]]
            push!(cargo_scenarios, CargoCollection(cargoes))
        elseif i == sc[2] # contract
            continue
        else
            # copy old scenario
            push!(cargo_scenarios, cargoC.items[i])
            probability[i] = pr[i]
        end
    end
    # remove extra element in probability
    new_probability = filter(x -> x!=0, probability)
    return CargoCollectionScenarios(cargo_scenarios), new_probability
end

# Naive scenario reduction. Comparing each cargo with it self in other scenarios
function Find_one_reduction_naive(CargoC::CargoCollectionScenarios, sc_old)
    n_cargo = length(CargoC.items[1])
    model = Model(Gurobi.Optimizer)
    set_silent(model) # removes terminal output
    #set_time_limit_sec(model, time_limit)
	@variable(model, x[1:sc_old,1:sc_old], Bin) # Binary variable for scenario reduction
    # Can only be contracted with one other scenario
    @constraint(model, [i = 1:sc_old], sum(x[i,j] for j ∈ 1:sc_old) <= 1)
    @constraint(model, [j = 1:sc_old], sum(x[i,j] for i ∈ 1:sc_old) <= 1)
    # Cannot be contracted with itself
    @constraint(model, [i = 1:sc_old], x[i,i] == 0) 
    # Contract two scenarios into 1
    @constraint(model, sum(x[i,j] for i ∈ 1:(sc_old-1), j ∈ (i+1):sc_old) == 1)
    # Cost of contraction of scenarios
    @expression(model, cost[i = 1:(sc_old-1),j = (i+1):sc_old],
        sum(abs.(getfield.(CargoC.items[i],:weight) .- getfield.(CargoC.items[j],:weight)))
        )
    # Objective function
    @objective(model, Min,
        sum(x[i,j]*cost[i,j] for i ∈ 1:(sc_old-1), j ∈ (i+1):sc_old)
        )
    optimize!(model)
    #println(cost)
    #termination_status(model)
    if termination_status(model)== MOI.OPTIMAL
        return value.(x), objective_value(model)
        
    else
        return 0,0
    end
end

#################################################

# Find two scenarios to contract. Compare scenarios and cargotype with eachother
# Heuristic
# Assume cargo of same type are interchangeable - not the case if considering hazardous cargo
function scenario_reduction_heuristic(CargoC::CargoCollectionScenarios, probability, sc, timelimit)
    n_types = [length(filter(x -> x.cargo_type_id == i, CargoC.items[1])) for i in 1:4] 
    no_scenarios = length(CargoC.items)
    it_timelimit = timelimit/(no_scenarios-sc)
    while no_scenarios > sc
        # Find two scenarios to contract
        # Find cost of moving cargo from scenario i to scenario id
        assignment, assignment_cost, scenarios_contracted = Find_one_reduction_heuristic(CargoC, it_timelimit)
        # contract the two scenarios
        CargoC, probability = contract_scenarios_heuristic(assignment, scenarios_contracted, CargoC, probability)
        #println("probability: ", probability)
        no_scenarios -= 1
    end
    return CargoC, probability
end

function Find_one_reduction_heuristic(CargoC::CargoCollectionScenarios, timelimit)
    sc = length(CargoC.items)
    #println("sc: ", sc)
    # randomly choose two scenarios to contract
    sc1 = rand(1:sc)
    sc2 = rand(1:sc)
    while sc1 == sc2
        sc2 = rand(1:sc)
    end
    #println(sc1, ", ", sc2)
    start = time_ns()
    assignment_best, reduction_cost_best = find_assignment(CargoC, sc1, sc2)
    scenarios_contracted = [sc1, sc2]
    elapsed = (time_ns() - start) / 1e9
    it = 1
    while elapsed < timelimit
        sc1 = rand(1:sc)
        sc2 = rand(1:sc)
        while sc1 == sc2
            sc2 = rand(1:sc)
        end
        # Find cost 
        assignment_temp, reduction_cost_temp = find_assignment(CargoC, sc1, sc2)
        if reduction_cost_temp < reduction_cost_best
            assignment_best = copy(assignment_temp)
            reduction_cost_best = reduction_cost_temp
            scenarios_contracted = [sc1, sc2]
        end
        elapsed = (time_ns() - start) / 1e9
        it += 1
    end
    #println("Did: ", it, " iterations to find one contraction")
    return assignment_best, reduction_cost_best, scenarios_contracted
end

function find_assignment(CargoC::CargoCollectionScenarios, sc1, sc2)
    # Build cost matrix
    n = length(CargoC.items[1])
    cost_matrix = Matrix{Union{Float64, Missing}}(undef, n, n)
    for i in 1:(n-1)
        for j in (i+1):n
            if (CargoC.items[sc1].items[i].cargo_type_id == CargoC.items[sc2].items[j].cargo_type_id)
                cost_matrix[i,j] = abs(CargoC.items[sc1].items[i].weight - CargoC.items[sc2].items[j].weight)
                cost_matrix[j,i] = cost_matrix[i,j] # symmetric matrix
            else
                cost_matrix[i,j] = missing
                cost_matrix[j,i] = missing
            end
        end
    end
    # Solve the assignment problem using Hungarian algorithm
    assignment, cost = hungarian(cost_matrix)
    return assignment, cost
end
# Contract the two scenarios
function contract_scenarios_heuristic(assignment, scenarios, 
    CargoC::CargoCollectionScenarios, probability)
    n = length(CargoC.items[1])
    sc = length(CargoC.items)
    new_prob = zeros(sc)
    cargo_scenarios = Vector{CargoCollection}()
    cargoes = Vector{Cargo}()
    for i in 1:sc
        if i == scenarios[1]
            # contract scenarios
            for j in 1:n
                # create new cargo
                current_cargo = CargoC.items[i].items[j]
                new_cargo = Cargo(
                    id=current_cargo.id,
                    cargo_type_id=current_cargo.cargo_type_id,  # Use the id from type_info
                    weight= (current_cargo.weight + CargoC.items[scenarios[2]].items[assignment[j]].weight)/2,
                    loading_port=current_cargo.loading_port,
                    discharge_port=current_cargo.discharge_port,
                    priority=current_cargo.priority,
                    requires_lashing=current_cargo.requires_lashing,
                    requires_ventilation=current_cargo.requires_ventilation,
                    hazardous=current_cargo.hazardous,
                    refers=current_cargo.refers
                )
                push!(cargoes, new_cargo)
            end
            new_prob[i] = probability[i] + probability[scenarios[2]]
            push!(cargo_scenarios, CargoCollection(cargoes))
        elseif i == scenarios[2]
            continue
        else
            # copy old scenario
            push!(cargo_scenarios, CargoC.items[i])
            new_prob[i] = probability[i]
        end
    end
    # remove extra element in probability
    new_probability = filter(x -> x!=0, new_prob)
    return CargoCollectionScenarios(cargo_scenarios), new_probability
end 

################################
# Uses clustering and K-means to reduce scenarios
function scenario_reduction_clustering(CargoC::CargoCollectionScenarios, probability, sc, timelimit)
    # Find Cost matrix
    sc_old = length(CargoC.items)
    cost = zeros(sc, sc)
    for i in 1:sc
        for j in 1:sc
            cost[i,j] = sum(abs.(getfield.(CargoC.items[i],:weight) .- getfield.(CargoC.items[j],:weight)))
        end
    end
    # Perform hierarchical clustering
    result = hclust(D, linkage = :average)  # or :single, :complete, :ward, etc.
    # Cut the dendrogram into 10 clusters
    labels = cutree(result, k = sc)
    if sum([count(x->x == i, labels) for i in 1:10]) != sc_old
        throw("Clustering did not produce the expected number of scenarios")
    end
    # Combine scenarios based on clustering labels
    for i in 1:length(labels)
        

    end



end
