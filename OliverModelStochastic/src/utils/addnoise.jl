# Script to add noise white noise to test instances.

# Snak med Oliver om det giver mening at gøre

# Adds white noise with mean 0 and standard deviation noise_level to the cargo weight
function add_white_noise_to_test_instance(problem::StowageProblem, noise_level::Float64)
    cargo = problem.cargo
    n = length(cargo)
    min_weight = [5.0 1.5 12.4 39.7] # From CARGO_TYPES in cargo.jl
    # Generate random noise for each cargo weight
    noise = randn(n) * noise_level
    # Add noise to the cargo weights
    new_cargoes = Vector{Cargo}()
    for i in 1:n
        current_cargo = cargo[i]
        new_cargo_weight = current_cargo.weight + noise[i]
        # Check if cargo is too light
        if new_cargo_weight < min_weight[current_cargo.cargo_type_id]
            new_cargo_weight = min_weight[current_cargo.cargo_type_id]
        end
        push!!(new_cargoes, Cargo(
            id = current_cargo.id,
            cargo_type_id = current_cargo.cargo_type_id,
            weight = new_cargo_weight,
            loading_port = current_cargo.loading_port,
            discharge_port = current_cargo.discharge_port,
            priority = current_cargo.priority,
            requires_lashing = current_cargo.requires_lashing,
            requires_ventilation = current_cargo.requires_ventilation,
            hazardous = current_cargo.hazardous,
            refers = current_cargo.refers))
    end
    new_cargocollection = CargoCollection(new_cargoes)
    # Return new problem
    new_problem = StowageProblem(
        vessel = problem.vessel,
        slots = problem.slots,
        cargo = new_cargocollection,    
        # Problem metadata
        name = problem.name,
        timestamp = problem.timestamp)
    return new_problem
end
