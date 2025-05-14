# Scenario reduction
function scenario_reduced(problem::StochasticStowageProblem,sc::Int64)
    cargoC = problem.cargo
    pr = problem.probability
    cargoc, probability = scenario_reduction(cargoC,pr, sc)
    # create new problem
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
function scenario_reduction(cargoC::CargoCollectionScenarios,probability::Vector{Float64}, sc::Int64)
    sc_old = length(cargoC.items)
    if sc_old <= sc
        throw("No reduction possible")
        return cargoC, probability
    end
    scenarios = sc_old
    while scenarios != sc
        x , obj_val = Find_one_reduction_naive(cargoC, scenarios)
        if x == 0
            println("No reduction possible")
            return cargoC
        end
        # Contract scenarios
        println("Reducing")
        println("")
        idx = findall(t -> t==1, x)
        sc1 = idx[1][1]
        sc2 = idx[1][2]
        cargoC, probability = contract_scenarios_naive(cargoC, probability, sc1, sc2, scenarios)
        println("Number of scenarios: ", length(cargoC.items))
        scenarios -= 1
    end
    return cargoC, probability
end

function contract_scenarios_naive(cargoC::CargoCollectionScenarios, pr::Vector{Float64},
    sc1::Int64, sc2::Int64, sc_old::Int64)
    # Create a new cargo collection with the reduced number of scenarios
    cargo_scenarios = Vector{CargoCollection}()
    cargoes = Vector{Cargo}()
    probability = zeros(sc_old) # array with probabilities
    for i in 1:sc_old
        if i == sc1 # contract
            for j in 1:length(cargoC.items[i].items)
                # create new cargo
                current_cargo = cargoC.items[i].items[j]
                new_cargo = Cargo(
                    id=current_cargo.id,
                    cargo_type_id=current_cargo.cargo_type_id,  # Use the id from type_info
                    weight=current_cargo.weight + cargoC.items[sc2].items[j].weight,
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
            probability[i] = pr[sc1] + pr[sc2]
            push!(cargo_scenarios, CargoCollection(cargoes))
        elseif i == sc2 # contract
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
        sum(abs(CargoC.items[i].items[c].weight - CargoC.items[j].items[c].weight) for c ∈ 1:n_cargo)
        )
    # Objective function
    @objective(model, Min,
        sum(x[i,j]*cost[i,j] for i ∈ 1:(sc_old-1), j ∈ (i+1):sc_old)
        )
    optimize!(model)
    #termination_status(model)
    if termination_status(model)== MOI.OPTIMAL
        return value.(x), objective_value(model)
        
    else
        return 0,0
    end
end

# NB! Not done
# Find two scenarios to contract. Compare scenarios and cargotype with eachother
# Assume cargo of same type are interchangeable - not the case if considering hazardous cargo
function Find_one_reduction(CargoC::CargoCollectionScenarios, sc_old)
    n_cargo = length(CargoC.items[1])

    # Find cost matrix and contraction for each index
    # TODO 
    contract = zeros(sc_old,sc_old)
    cost = zeros(sc_old,sc_old)
    for i in 1:(sc_old-1), j in i+1:sc_old
        # TODO: Solve function
        contract[i,j] = 0
        cost[i,j] = 0
    end 
    # Find which should be contracted
    model = Model(Gurobi.Optimizer)
    set_silent(model) # removes terminal output
    #set_time_limit_sec(model, time_limit)
    @variable(model, x[1:sc_old,1:sc_old], Bin) # Binary variable for scenario reduction
    @constraint(model, [i = 1:sc_old], sum(x[i,j] for j ∈ 1:sc_old) <= 1)
    @constraint(model, [j = 1:sc_old], sum(x[i,j] for i ∈ 1:sc_old) <= 1)
    # Cannot be contracted with itself
    @constraint(model, [i = 1:sc_old], x[i,i] == 0) 
    # Contract two scenarios into 1
    @constraint(model, sum(x[i,j] for i ∈ 1:(sc_old-1), j ∈ (i+1):sc_old) == 1)
    # Cost of contraction of scenarios
    @objective(model, Min, sum(cost[i,j]*x[i,j] for i ∈ 1:(sc_old-1), j ∈ (i+1):sc_old))
    optimize!(model)

    if termination_status(model)== MOI.OPTIMAL
        return value.(x), contract
    else
        return 0,0
    end
end 
# Find the cheapest way to combine two scenarios, where cargo only 
# has to match cargo of similar type
# TODO: Should be a heuristic
function contraction_cost(cargo1, cargo2)
    println("TODO")
end


