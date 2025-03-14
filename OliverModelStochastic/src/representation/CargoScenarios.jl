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
function generate_simple_cargo_scenarios(cargoC::CargoCollection,s::Int64, ids)
        cargo_scenarios = Vector{CargoCollection}() # array with CargoCollections 
        probability = fill(1/s,s) # equal probability for each scenario
        for i in 1:s
                push!(cargo_scenarios,random_weight_CarcoCollection(cargoC,ids))
        end
        return CargoCollectionScenarios(cargo_scenarios), probability # return CargoCollectionScenarios
end

###########################

# Finds mean and variance of each cargo type
function mean_var_car(problem::StowageProblem)
        cargo1 = filter(x -> x.id == 1,problem.cargo.items) # find cargo fo type 1
        cargo2 = filter(x -> x.id == 2,problem.cargo.items) # find cargo fo type 2
        cargo3 = filter(x -> x.id == 3,problem.cargo.items) # find cargo fo type 3
        cargo4 = filter(x -> x.id == 4,problem.cargo.items) # find cargo fo type 4
        weight1 = [cargo.weight for cargo in cargo1] # get weights
        weight2 = [cargo.weight for cargo in cargo2]
        weight3 = [cargo.weight for cargo in cargo3]
        weight4 = [cargo.weight for cargo in cargo4]
        means = [mean(weight1),mean(weight2),mean(weight3),mean(weight4)] # calculate mean
        variances = [var(weight1),var(weight2),var(weight3),var(weight4)] # calculate variance
        return means,variances
end
# Monto Carlo sampling for cargo weights. Weight for each cargo is assumed normal distributed
function Monto_Carlo_sampling()

end

# Needs to return CargoCollectionScenarios() and probability


