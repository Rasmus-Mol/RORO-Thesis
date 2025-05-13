@kwdef struct StowageProblem
    # Core components
    vessel::Vessel
    slots::SlotCollection
    cargo::CargoCollection    
    # Problem metadata
    name::String
    timestamp::DateTime = now()
end

# Generate StochasticStowageProblem
@kwdef struct StochasticStowageProblem
    # Core components
    vessel::Vessel
    slots::SlotCollection
    cargo::CargoCollectionScenarios   
    unknown_weights::Vector{Int64}
    known_weights::Vector{Int64}
    scenarios::Int64
    probability::Vector{Float64}
    # Problem metadata
    name::String
    timestamp::DateTime = now()
end

function load_data(vessel_name::String, instance_name::String, instance_type::String)
    # Load vessel data
    vessel = Vessel(vessel_name)
    vessel = simplify_vessel(vessel, target_points=50)

    
    # Load slots from data directory
    slots = load_all_slots(joinpath(@__DIR__, "../../data/slots", vessel_name), vessel_name)
    
    # Load cargo and metadata from default instance
    instance_path = joinpath(@__DIR__, "../../instances/", instance_type, vessel_name, "$instance_name.json")
    cargo, metadata = load_instance(instance_path)
    
    return StowageProblem(
        vessel = vessel,
        slots = slots,
        cargo = cargo,
        name = metadata.name
    )
end

# Generates the stochastic problem
# Generate CargoCollectionScenarios. s is number of scenarios, n is number of random cargos with weight changed
# generate_method is a function that generates scenarios, default is generate_simple_cargo_scenarios
function create_stochastic_problem(problem::StowageProblem, s::Int64,n::Int64, ids = [], generate_method::Function = generate_simple_cargo_scenarios)
    ids = ids == [] ? random_ids(problem.cargo,n) : ids # if ids are not given, generate random ids
    cargo_scenarios, probability = generate_method(problem, s, ids)
    return StochasticStowageProblem(
        vessel = problem.vessel,
        slots = problem.slots,
        cargo = cargo_scenarios,
        unknown_weights = sort(ids),
        known_weights = sort(setdiff([i for i in 1:length(problem.cargo.items)],ids)),
        scenarios = s,
        probability = probability,
        name = problem.name
    )
end


# Create expected value problem from stochastic problem
function expected_value_problem(problem::StochasticStowageProblem)
    n_cargo = length(problem.cargo.items[1])
    cargos = Vector{Cargo}()
    # Find mean weight for each cargo
    for i in 1:n_cargo
        weight = 0
        for j in 1:problem.scenarios
            weight += problem.cargo.items[j].items[i].weight
        end
        push!(cargos, Cargo(
            id = problem.cargo.items[1].items[i].id,
            cargo_type_id = problem.cargo.items[1].items[i].cargo_type_id,
            weight = weight/problem.scenarios,
            loading_port = problem.cargo.items[1].items[i].loading_port,
            discharge_port = problem.cargo.items[1].items[i].discharge_port,
            priority = problem.cargo.items[1].items[i].priority,
            requires_lashing = problem.cargo.items[1].items[i].requires_lashing,
            requires_ventilation = problem.cargo.items[1].items[i].requires_ventilation,
            hazardous = problem.cargo.items[1].items[i].hazardous,
            refers = problem.cargo.items[1].items[i].refers
        ))
    end
    return StowageProblem(
        vessel = problem.vessel,
        slots = problem.slots,
        cargo = CargoCollection(cargos),
        name = problem.name
    )
end