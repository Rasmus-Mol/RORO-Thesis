# Scenario reduction
function scenario_reduction(cargoC::CargoCollectionScenarios,probability::Vector{Float64}, sc::Int64)
    sc_old = length(cargoC.items)
    x = asd
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
        sum(abs(cargoC.items[i].items[c].weight - cargoC.items[j].items[c].weight) for c ∈ 1:n_cargo)
        )
    # Objective function
    @objective(model, Min,
        sum(x[i,j]*cost[i,j] for i ∈ 1:(sc_old-1), j ∈ (i+1):sc_old)
        )
    optimize!(model)
    return value.(x), objective_value(model)
end

# Find two scenarios to contract. Compare scenarios and cargotype with eachother
# Assume cargo of same type are interchangeable - not the case if considering hazardous cargo
function Find_one_reduction(CargoC::CargoCollectionScenarios, sc_old)
    n_cargo = length(CargoC.items[1])
    model = Model(Gurobi.Optimizer)
    set_silent(model) # removes terminal output
    #set_time_limit_sec(model, time_limit)
    # Number of each cargo-type
    n_types = [length(filter(x -> x.cargo_type_id == i, cargoC.items[1])) for i in 1:4]
    # Binary variables for each cargo type
    # x1[i,j,n1,n2] = 1 if cargo n1 in scenario i is matched with cargo n2 in scenario j
	@variable(model, x1[1:sc_old,1:sc_old,1:n_types[1],1:n_types[1]], Bin) # Binary variable for type 1
    @variable(model, x2[1:sc_old,1:sc_old,1:n_types[2],1:n_types[1]], Bin) # Binary variable for type 2
    @variable(model, x3[1:sc_old,1:sc_old,1:n_types[3],1:n_types[1]], Bin) # Binary variable for type 3
    @variable(model, x4[1:sc_old,1:sc_old,1:n_types[4],1:n_types[1]], Bin) # Binary variable for type 4

    # cargo have to be matched to a cargo in different scenario
    @constraint(model, [c1=1:n_types[1]], sum(x1[i,j,c1,c2] for c2=1:n_types[1], i = 1:(sc_old-1), j=(i+1):sc_old) == 1)
    @constraint(model, [c1=1:n_types[2]], sum(x2[i,j,c1,c2] for c2=1:n_types[1], i = 1:(sc_old-1), j=(i+1):sc_old) == 1)
    @constraint(model, [c1=1:n_types[3]], sum(x3[i,j,c1,c2] for c2=1:n_types[1], i=1:(sc_old-1), j=(i+1):sc_old) == 1)
    @constraint(model, [c1=1:n_types[4]], sum(x4[i,j,c1,c2] for c2=1:n_types[1], i=1:(sc_old-1), j=(i+1):sc_old) == 1)
    
    # only match two scenarios - probably not needed
    #=@constraint(model, [i = 1:sc_old], sum(x1[i,j,n] for j ∈ 1:sc_old, n ∈ 1:n_types[1]) <= 1)
    @constraint(model, [j = 1:sc_old], sum(x1[i,j,n] for i ∈ 1:sc_old, n ∈ 1:n_types[1]) <= 1)
    @constraint(model, [i = 1:sc_old], sum(x2[i,j,n] for j ∈ 1:sc_old, n ∈ 1:n_types[2]) <= 1)
    @constraint(model, [j = 1:sc_old], sum(x2[i,j,n] for i ∈ 1:sc_old, n ∈ 1:n_types[2]) <= 1)
    @constraint(model, [i = 1:sc_old], sum(x3[i,j,n] for j ∈ 1:sc_old, n ∈ 1:n_types[3]) <= 1)
    @constraint(model, [j = 1:sc_old], sum(x3[i,j,n] for i ∈ 1:sc_old, n ∈ 1:n_types[3]) <= 1)
    @constraint(model, [i = 1:sc_old], sum(x4[i,j,n] for j ∈ 1:sc_old, n ∈ 1:n_types[4]) <= 1)
    @constraint(model, [j = 1:sc_old], sum(x4[i,j,n] for i ∈ 1:sc_old, n ∈ 1:n_types[4]) <= 1)
    =#

    # Cannot contract scenario with itself
    @constraint(model, [i = 1:sc_old,c1=1:n_types[1],c2=1:n_types[1]], x1[i,i,c1,c2] == 0)
    @constraint(model, [i = 1:sc_old,c1=1:n_types[2],c2=1:n_types[1]], x2[i,i,c1,c2] == 0)
    @constraint(model, [i = 1:sc_old,c1=1:n_types[3],c2=1:n_types[1]], x3[i,i,c1,c2] == 0)
    @constraint(model, [i = 1:sc_old,c1=1:n_types[4],c2=1:n_types[1]], x4[i,i,c1,c2] == 0)

    # Cost for each type
    @expression(model, cost1[i = 1:sc_old,j = 1:sc_old,n = 1:n_types[1]],
        sum(abs(cargoC.items[i].items[n].weight - cargoC.items[j].items[c.id].weight) for c ∈ filter(x -> x.cargo_type_id == 1, cargoC.items[1]))
        )
    @expression(mode, cost2[i = 1:sc_old,j = 1:sc_old,n = 1:n_types[2]],
        sum(abs(cargoC.items[i].items[n].weight - cargoC.items[j].items[c.id].weight) for c ∈ filter(x -> x.cargo_type_id == 2, cargoC.items[1]))
        )
    @expression(model, cost3[i = 1:sc_old,j = 1:sc_old,n = 1:n_types[3]],
        sum(abs(cargoC.items[i].items[n].weight - cargoC.items[j].items[c.id].weight) for c ∈ filter(x -> x.cargo_type_id == 3, cargoC.items[1]))
        )
    @expression(model, cost4[i = 1:sc_old,j = 1:sc_old,n = 1:n_types[4]],
        sum(abs(cargoC.items[i].items[n].weight - cargoC.items[j].items[c.id].weight) for c ∈ filter(x -> x.cargo_type_id == 4, cargoC.items[1]))
        )
    # Objective function
    @objective(model,Min,
    sum(x1[i,j]*cost[i,j] for i ∈ 1:(sc_old-1), j ∈ (i+1):sc_old)
    )
end

#sto_pro = create_stochastic_problem(problem, scenarios, n_cargo_unknownweight,[])
cargoC = sto_pro.cargo
x, obj = Find_one_reduction_naive(cargoC, 20)
x
indices = findall(x .== 1)

# test model
sc_old = 20
model = Model(Gurobi.Optimizer)
n_types = [length(filter(x -> x.cargo_type_id == i, cargoC.items[1])) for i in 1:4]
@variable(model, x1[1:sc_old,1:sc_old,1:n_types[1],1:n_types[1]], Bin) # Binary variable for type 1
    @variable(model, x2[1:sc_old,1:sc_old,1:n_types[2],1:n_types[1]], Bin) # Binary variable for type 2
    @variable(model, x3[1:sc_old,1:sc_old,1:n_types[3],1:n_types[1]], Bin) # Binary variable for type 3
    @variable(model, x4[1:sc_old,1:sc_old,1:n_types[4],1:n_types[1]], Bin) # Binary variable for type 4

test = x1[]


