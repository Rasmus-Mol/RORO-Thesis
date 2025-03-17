# To generate scenarios and put them into struct for problem
using StatsBase
# Randomly select n number of cargos to have random weight
function random_ids(cargo::CargoCollection,n)
        ids = cargo.items.id
        return sample(ids,n,replace = false)
end

# Random weight depending on type. Range is from cargo.jl
function random_weight(type)
        Rweight = 0
        if type == 2 # Car
                Rweight = rand(1.5:0.1:3.0)
        elseif type == 1 # Truck
                Rweight = rand(5.0:0.5:40.0)
        elseif type == 3 # Machine
                volume = 4.6 * 4.5 * 3.0
                Rweight = volume * rand(0.2:0.01:0.5)
        else # Secu
                volume = 13.8 * 3.6 * 4.0
                Rweight = volume * rand(0.2:0.01:0.5)
        end
        return Rweight
end
# Generates new CargoCollection with random weight for given ids
function random_weight_CarcoCollection(cargoC::CargoCollection,ids) 
        new_cargoC = Vector{Cargo}() # new CargoCollection
        n_items = length(cargoC.items) # number of cargo 
        for i in 1:n_items
                current_cargo = cargoC.items[i]
                if current_cargo.id in ids # generate new cargo with random weight
                        push!(new_cargoC,Cargo(id = current_cargo.id,
                        cargo_type_id = current_cargo.cargo_type_id,
                        weight = random_weight(current_cargo.cargo_type_id),
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
        return CargoCollection(new_cargoC) # return new CargoCollection with some weights changed
end
# Generate CargoCollectionScenarios. s is number of scenarios
function generate_simple_cargo_scenarios(problem::StowageProblem,s::Int64, ids)
        cargo_scenarios = Vector{CargoCollection}() # array with CargoCollections 
        probability = fill(1/s,s) # equal probability for each scenario
        cargoC = problem.cargo
        for i in 1:s
                push!(cargo_scenarios,random_weight_CarcoCollection(cargoC,ids))
        end
        return CargoCollectionScenarios(cargo_scenarios), probability # return CargoCollectionScenarios
end

###########################

# Finds mean and variance of each cargo type
function mean_var_car(problem::StowageProblem)
        # Get weights for each type car
        weight1 = [cargo.weight for cargo in filter(x -> x.cargo_type_id == 1,problem.cargo.items)] 
        weight2 = [cargo.weight for cargo in filter(x -> x.cargo_type_id == 2,problem.cargo.items)] 
        weight3 = [cargo.weight for cargo in filter(x -> x.cargo_type_id == 3,problem.cargo.items)] 
        weight4 = [cargo.weight for cargo in filter(x -> x.cargo_type_id == 4,problem.cargo.items)] 
        weights = [weight1,weight2,weight3,weight4]
        means = []
        variances = []
        for w in weights  # 4 types of cargo 
                if length(w)>0
                        push!(means,mean(w)) # calculate mean
                        push!(variances,var(w)) # calculate variance
                else
                        push!(means,-1)
                        push!(variances,-1)
                end
        end
        return means,variances
end
# Monto Carlo sampling for cargo weights. Weight for each cargo type is assumed normal distributed
function Monto_Carlo_sampling(problem::StowageProblem,sc::Int64,ids)
        means,variances = mean_var_car(problem) # get means and variances
        # Normal distributions for each cargo type
        distributions = []
        for i in 1:4
                if means[i] != -1 # -1 indicates that there is no cargo of that type
                        push!(distributions,Normal(means[i],sqrt(variances[i])))
                else
                        push!(distributions,-1)
                end
        end
        n_items = length(problem.cargo.items) # number of cargos
        #sample_weight = 0
        cargo_scenarios = Vector{CargoCollection}() # array with CargoCollections

        for i in 1:sc # number of scenarios to generate
                new_cargoC = Vector{Cargo}() # new CargoCollection
                for j in 1:n_items
                        current_cargo = problem.cargo.items[j]
                        if j in ids # generate random weight
                                sample_weight = rand(distributions[problem.cargo.items[j].cargo_type_id]) # generate random weight
                                push!(new_cargoC,Cargo(id = current_cargo.id,
                                cargo_type_id = current_cargo.cargo_type_id,
                                weight = sample_weight,
                                loading_port = current_cargo.loading_port,
                                discharge_port = current_cargo.discharge_port,
                                priority = current_cargo.priority,
                                requires_lashing = current_cargo.requires_lashing,
                                requires_ventilation = current_cargo.requires_ventilation,
                                hazardous = current_cargo.hazardous,refers = current_cargo.refers))
                        else
                                push!(new_cargoC,current_cargo) # copy cargo
                        end
                end
                push!(cargo_scenarios,CargoCollection(new_cargoC)) # add new CargoCollection to array
        end
        return CargoCollectionScenarios(cargo_scenarios), fill(1/sc,sc) # return CargoCollectionScenarios
end

# Needs to return CargoCollectionScenarios() and probability


