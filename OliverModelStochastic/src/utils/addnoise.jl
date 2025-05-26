# Script to add noise white noise to test instances.

# Adds white noise with mean 0 and standard deviation noise_level to the cargo weight
# Default noise is [truck, car, machine, secu] = []
function add_white_noise_to_test_instance(problem::StowageProblem, noise_level::Vector{Float64} = [2.5,0.1,1.2,4])
    cargo = problem.cargo
    n = length(cargo)
    n_types = [length(filter(c -> c.cargo_type_id == i, cargo)) for i in 1:4]
    noise_types = ones(n_types[1]) * noise_level[1]
    for i in 2:4
        noise_types = vcat(noise_types, ones(n_types[i]) * noise_level[i])
    end
    min_weight = [5.0 1.5 12.4 39.7] # From CARGO_TYPES in cargo.jl
    max_weight = [40, 3, 30.375, 99.36] # From CARGO_TYPES in cargo.jl
    # Generate random noise for each cargo weight
    noise = randn(n) .* noise_types
    idx = [1,1,1,1] # keep track of how many of each has been used
    # Add noise to the cargo weights
    new_cargoes = Vector{Cargo}()
    for i in 1:n
        current_cargo = cargo[i]
        if current_cargo.cargo_type_id == 1
            new_cargo_weight = current_cargo.weight + noise[idx[1]]
            idx[1] += 1
        elseif current_cargo.cargo_type_id == 2
            new_cargo_weight = current_cargo.weight + noise[n_types[1] + idx[2]]
            idx[2] += 1
        elseif current_cargo.cargo_type_id == 3
            new_cargo_weight = current_cargo.weight + noise[n_types[1]+n_types[2] + idx[3]]
            idx[3] += 1
        elseif current_cargo.cargo_type_id == 4
            new_cargo_weight = current_cargo.weight + noise[n_types[1]+n_types[2]+n_types[3] + idx[4]]
            idx[4] += 1
        end
        # Check if cargo is too light
        if new_cargo_weight < min_weight[current_cargo.cargo_type_id]
            new_cargo_weight = min_weight[current_cargo.cargo_type_id]
        elseif new_cargo_weight > max_weight[current_cargo.cargo_type_id]
            new_cargo_weight = max_weight[current_cargo.cargo_type_id]
        end
        push!(new_cargoes, Cargo(
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
