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

# Adds biased noise, ie. so it corresponds to the historic trend
function add_biased_noise_to_test_instance(problem::StowageProblem)
    CargoC = problem.cargo
    ids = [i for i in 1:length(CargoC.items)]
    cur_path = @__DIR__
    file_path = joinpath(cur_path,"..","..","data","CargoWeights.csv")
    #file_path = joinpath(cur_path,"data","CargoWeights.csv")
    # load historic data 
    df = load_Weight_Variance_data(file_path)
    df_secu, df_trailer = weight_difference_info(df,false)
    # Seperate data into quantiles
    n_quantiles = 4
    df_secu, q_arr_secu = seperate_data_into_quantiles(df_secu,n_quantiles,false)
    df_trailer, q_arr_trailer = seperate_data_into_quantiles(df_trailer,n_quantiles,false)
    q_arr_secu = q_arr_secu ./1000 # convert to tons
    q_arr_trailer = q_arr_trailer ./1000 # convert to tons
    secu_var = [collect(filter(x -> x.QuantileNumber == i, df_secu).Variance) for i in 1:n_quantiles]
    trailer_var = [collect(filter(x -> x.QuantileNumber == i, df_trailer).Variance) for i in 1:n_quantiles]
    new_cargoc = biased_noise_quantile_cargocollection(CargoC,ids,q_arr_secu,q_arr_trailer,secu_var,trailer_var)
    # Return new problem
    return new_problem = StowageProblem(
        vessel = problem.vessel,
        slots = problem.slots,
        cargo = new_cargoc,
        # Problem metadata
        name = problem.name,
        timestamp = problem.timestamp)
end
# Generate biased noise to 1 scenario
function biased_noise_quantile_cargocollection(cargo_C,ids,q_arr_secu,q_arr_trailer,secu_var,trailer_var)
    new_cargoC = Vector{Cargo}() # new CargoCollection
    n_items = length(cargo_C.items) # number of cargo 
    for i in 1:n_items
            current_cargo = cargo_C.items[i]
            if current_cargo.id in ids # generate new cargo with new weight
                    push!(new_cargoC,Cargo(id = current_cargo.id,
                    cargo_type_id = current_cargo.cargo_type_id,
                    weight =  biased_noise_bookedweight_quantile_cargovariance(current_cargo,q_arr_secu,q_arr_trailer,secu_var,trailer_var),
                    loading_port = current_cargo.loading_port,
                    discharge_port = current_cargo.discharge_port,
                    priority = current_cargo.priority,
                    requires_lashing = current_cargo.requires_lashing,
                    requires_ventilation = current_cargo.requires_ventilation,
                    hazardous = current_cargo.hazardous,refers = current_cargo.refers))
            else # copy cargo
                    push!(new_cargoC,current_cargo)
            end
    end
    return CargoCollection(new_cargoC)
end
# Generate biased noise for 1 cargo
function biased_noise_bookedweight_quantile_cargovariance(cargo::Cargo,
    q_arr_secu,q_arr_trailer,secu_var,trailer_var)
    c_type = cargo.cargo_type_id
    if c_type == 1 # Truck/Trailer
        for i in 1:(length(trailer_var)-1)
            if cargo.weight <= q_arr_trailer[i]
                # uniform sampling from variance from corresponding quantile
                new_weight = cargo.weight+rand(trailer_var[i])/1000 # convert to tons
                if new_weight < 5.0
                    return 5.0
                end
                return new_weight
            end
        end
        # Belongs to last quantile
        new_weight = cargo.weight+rand(trailer_var[end])/1000
        if new_weight < 5.0
            return 5.0
        end
        return new_weight
    elseif c_type == 2 # Car - no historical data
            return rand(1.5:0.1:3.0)
    elseif c_type == 3 # Machine - no historical data
            volume = 4.6 * 4.5 * 3.0
            return volume * rand(0.2:0.01:0.5)
    else # Secu
        # Checking which quantile the cargo weight is in
        for i in 1:(length(secu_var)-1)
            if cargo.weight <= q_arr_secu[i]
                # uniform sampling from variance from corresponding quantile
                new_weight =  cargo.weight+rand(secu_var[i])/1000 # convert to tons
                if new_weight < 13.8 * 3.6 * 4.0 * 0.2
                    return 13.8 * 3.6 * 4.0 * 0.2 # min weight
                end
                return new_weight
            end
        end
        # Belongs to last quantile
        new_weight = cargo.weight+rand(secu_var[end])/1000
        if new_weight < 13.8 * 3.6 * 4.0 * 0.2
            return 13.8 * 3.6 * 4.0 * 0.2 # min weight
        end
        return new_weight
    end
end
