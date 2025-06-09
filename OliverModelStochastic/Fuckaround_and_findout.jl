include("packages_and_files.jl")

# Do stuff 
test_instance_hol = Hollandia_test[5]
test_instance_fin = Finlandia_test[5]

Deterministic_problem_hol = load_data("hollandia",test_instance_hol,"hazardous")
#Deterministic_problem_fin = load_data("finlandia","no_cars_medium_100_haz_eq_0.1","hazardous")
Deterministic_problem_fin = load_data("finlandia",test_instance_fin,"hazardous")

Deterministic_problem_hol.vessel.decks
slots_deck1 = filter(x -> x.deck_id == 1, Deterministic_problem_hol.slots)
slots_deck1_truck = filter(x -> x.cargo_type_id == 1, slots_deck1)
slots_deck1_car = filter(x -> x.cargo_type_id == 2, slots_deck1)
slots_deck1_heavy = filter(x -> x.cargo_type_id == 3, slots_deck1)
slots_deck1_secu = filter(x -> x.cargo_type_id == 4, slots_deck1)
length(slots_deck1_truck)
length(slots_deck1_car)
length(slots_deck1_heavy)
length(slots_deck1_secu)

pro_hol = create_stochastic_problem(Deterministic_problem_hol, 10, length(Deterministic_problem_hol.cargo), []) 
pro_fin = create_stochastic_problem(Deterministic_problem_fin, 100, length(Deterministic_problem_fin.cargo), [])
pro_fin_1 = create_stochastic_problem(Deterministic_problem_fin, 10, length(Deterministic_problem_fin.cargo), [], Bootstrap_bookedweight_quantile)

pro_fin_2 = create_stochastic_problem_scenarioreduction(Deterministic_problem_fin, 10,
length(Deterministic_problem_fin.cargo), 100, [], Bootstrap_bookedweight_quantile, scenario_reduction_heuristic, 5*60)

pro_fin_3 = create_stochastic_problem_scenarioreduction(Deterministic_problem_fin, 10,
length(Deterministic_problem_fin.cargo), 50, [], Bootstrap_bookedweight_quantile, scenario_reduction_naive, 60)



CargoC = pro_fin.cargo
pro_fin_3.probability
sc_old = 100
cost = zeros(sc_old, sc_old)
for i in 1:(sc_old-1)
    for j in (i+1):sc_old
        cost[i,j] = sum(abs.(getfield.(CargoC.items[i],:weight) .- getfield.(CargoC.items[j],:weight)))
    end
end
cost2 = zeros(sc_old, sc_old)
for i in 1:sc_old
    for j in 1:sc_old
        cost2[i,j] = sum(abs.(getfield.(CargoC.items[i],:weight) .- getfield.(CargoC.items[j],:weight)))
    end
end

heatmap(cost)
heatmap(cost2)
using Clustering
temp = kmeans(cost, 10)
println("Cluster assignments: ", temp.assignments)
println("Centroids: ", temp.centers)
temp2 = kmeans(cost2, 10)
println("Cluster assignments: ", temp2.assignments)
println("Centroids: ", temp2.centers)

# Perform hierarchical clustering
result = hclust(cost2, linkage = :average)  # or :single, :complete, :ward, etc.
# Cut the dendrogram into 10 clusters
labels = cutree(result, k = 10)
println("Cluster assignments: ", labels)
println([count(x->x == i, labels) for i in 1:10])
sum([count(x->x == i, labels) for i in 1:10])



plot(cargoc.items.total_weight)
plot(pro_fin_3.cargo.items.total_weight)
plot!(pro_fin_2[1].cargo.items.total_weight)


plot(pro_fin_3.probability, label = "Naive", title = "initial prob = $(1/sc),\n initial scenarios = $(sc)")
plot!(pro_fin_2[1].probability, label = "Heuristic")



ids = random_ids(Deterministic_problem_fin.cargo,length(Deterministic_problem_fin.cargo))
sc = 200
cargo_scenarios, probability = Bootstrap_bookedweight_quantile(Deterministic_problem_fin, sc, random_ids(Deterministic_problem_fin.cargo,length(Deterministic_problem_fin.cargo)))
#generate_simple_cargo_scenarios(Deterministic_problem_fin, sc, random_ids(Deterministic_problem_fin.cargo,length(Deterministic_problem_fin.cargo)))
weight_id_1 = [cargo_scenarios.items[i].items[1].weight for i in 1:sc]

CargoC1, probability1 = scenario_reduction_heuristic(cargo_scenarios, probability, 20, 5*60)
weight_id_1_H = [CargoC1.items[i].items[1].weight for i in 1:5]

CargoC2, probability2 = scenario_reduction_naive(cargo_scenarios,probability, 20, 5*60)
weight_id_1_naive = [CargoC2.items[i].items[1].weight for i in 1:5]

println(probability2)
println(probability1)
plot(probability2, label = "Naive", title = "initial prob = $(1/sc),\n initial scenarios = $(sc)")
plot!(probability1, label = "Heuristic")


p1 = plot(Deterministic_problem_fin.cargo.items[1].weight*ones(sc), label = "Original weight")
scatter!(weight_id_1_H, label = "Heuristic")
scatter!(weight_id_1, label = "Original scenarios")

p2 = plot(Deterministic_problem_fin.cargo.items[1].weight*ones(sc), label = "Original weight")
scatter!(weight_id_1_naive, label = "Naive")
scatter!(weight_id_1, label = "Original scenarios")

cost1 = zeros(sc, sc)
cost2 = zeros(sc, sc)
for i in 1:(sc-1), j in (i+1):sc
    cost1[i,j] = sum(abs(cargo_scenarios.items[i].items[c].weight - cargo_scenarios.items[j].items[c].weight) for c in 1:length(cargo_scenarios.items[1].items))
    cost2[i,j] = sum(abs.(getfield.(cargo_scenarios.items[i],:weight) .- getfield.(cargo_scenarios.items[j],:weight)))
end
idx = findall(x -> x > 0.5, cost2)
vals = cost[idx]
min_val, relative_idx = findmin(vals)
idx[relative_idx]
cost2[97,115]

assignment, cost = find_assignment(cargo_scenarios, 97, 155)
n = length(cargo_scenarios.items[1])
    cost_matrix = Matrix{Union{Float64, Missing}}(undef, n, n)
    for i in 1:n
        for j in 1:n
            if (cargo_scenarios.items[97].items[i].cargo_type_id == cargo_scenarios.items[155].items[j].cargo_type_id)
                cost_matrix[i,j] = abs(cargo_scenarios.items[97].items[i].weight - cargo_scenarios.items[155].items[j].weight)
            else
                cost_matrix[i,j] = missing
            end
        end
    end
cost_matrix[1,5]
cost_matrix[5,1]
abs(cargo_scenarios.items[97].items[1].weight - cargo_scenarios.items[155].items[5].weight)
abs(cargo_scenarios.items[155].items[5].weight - cargo_scenarios.items[97].items[1].weight)
cargo_scenarios.items[97].items[5].weight
cargo_scenarios.items[155].items[1].weight
abs(10.5-35)
#######################
# Hollandia
model_hol = create_model(Deterministic_problem_hol)
set_silent(model_hol) # removes terminal output
set_time_limit_sec(model_hol, 60 * 15) # 5 minutes to solve model
optimize!(model_hol)
using MathOptInterface
const MOI = MathOptInterface
if termination_status(model_hol) in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
    println("Model is infeasible. Computing conflict...")
    MOI.compute_conflict!(JuMP.backend(model_hol))
end
for (func, set) in JuMP.list_of_constraint_types(model_hol)
    for constraint_ref in JuMP.all_constraints(model_hol, func, set)
        status = MOI.get(
            JuMP.backend(model_hol),
            MOI.ConstraintConflictStatus(),
            JuMP.constraint_object(constraint_ref)
        )
        if status == MOI.IN_CONFLICT
            println("Constraint in conflict: ", constraint_ref)
        end
    end
end

for v in JuMP.all_variables(model_hol)
    var_index = MOI.VariableIndex(v)
    lb_status = MOI.get(JuMP.backend(model_hol), MOI.VariableConflictStatus(), var_index, MOI.LOWER)
    ub_status = MOI.get(JuMP.backend(model_hol), MOI.VariableConflictStatus(), var_index, MOI.UPPER)

    if lb_status == MOI.IN_CONFLICT
        println("Lower bound of variable $(v) is in conflict")
    end
    if ub_status == MOI.IN_CONFLICT
        println("Upper bound of variable $(v) is in conflict")
    end
end





model_stochastic_hol = create_model_stochastic(pro_hol)
set_silent(model_stochastic_hol) # removes terminal output
set_time_limit_sec(model_stochastic_hol, 60 * 2) # 5 minutes to solve model
#set_time_limit_sec(model_stochastic, 60 * 15) # 10 minutes to solve model
optimize!(model_stochastic_hol)

sol_hol = extract_stochastic_solution(pro_hol,model_stochastic_hol)
folder1 = "Temp1_delete"
folder2 = "Temp2_delete"
write_solution_stochastic(sol_hol,folder2,"Stochastic_Solution",folder1)
cs_sol = sol_hol.cs
second_stage_m = second_stage_model(cs_sol, Deterministic_problem_hol)
set_silent(second_stage_m) # removes terminal output
set_time_limit_sec(second_stage_m, 60) 
optimize!(second_stage_m)
fitted_sol = get_solution_second_stage_stochastic(Deterministic_problem_hol, second_stage_m, sol_hol)
write_solution(fitted_sol,folder1,"Fitted_Solution",folder2)

#########################
# Finlandia 
model_fin = create_model(Deterministic_problem_fin)
set_silent(model_fin) # removes terminal output
set_time_limit_sec(model_fin, 60 * 5) # 5 minutes to solve model
optimize!(model_fin)
solution_det = extract_solution(Deterministic_problem_fin, model_fin)

# med ingen ændringer
# 89 sek med korrekt vcg slopes
# 41 sek med korrekt z_min og z_max
# med begge ændringer

model_stochastic_fin = create_model_stochastic(pro_fin)
set_silent(model_stochastic_fin) # removes terminal output
set_time_limit_sec(model_stochastic_fin, 60 * 5) # 5 minutes to solve model
#set_time_limit_sec(model_stochastic, 60 * 15) # 10 minutes to solve model
optimize!(model_stochastic_fin)
sol = extract_stochastic_solution(pro_fin,model_stochastic_fin)
cs_sol = sol.cs

# Solve second stage when we know unknown weights
second_stage_m = second_stage_model(cs_sol, problem_det)
set_silent(second_stage_m) # removes terminal output
set_time_limit_sec(second_stage_m, time_limit) 
optimize!(second_stage_m)

#=
for i in 1:length(vessel.decks)
    println("Deck $(i) weight limit: ", vessel.decks[i].weight_limit)
end

deck1 = filter(x -> x.deck_id == 1, slots)
slots1 = [filter(x->x.cargo_type_id == i, deck1) for i in 1:4]
println("Deck 1")
println("trucks: ", length(slots1[1]))
println("cars: ", length(slots1[2]))
println("Heavy machinery: ", length(slots1[3]))
println("Secu: ", length(slots1[4]))

deck2 = filter(x -> x.deck_id == 2, slots)
slots2 = [filter(x->x.cargo_type_id == i, deck2) for i in 1:4]
println("Deck 2")
println("trucks: ", length(slots2[1]))
println("cars: ", length(slots2[2]))
println("Heavy machinery: ", length(slots2[3]))
println("Secu: ", length(slots2[4]))

deck3 = filter(x -> x.deck_id == 3, slots)
slots3 = [filter(x->x.cargo_type_id == i, deck3) for i in 1:4]
println("Deck 3")
println("trucks: ", length(slots3[1]))
println("cars: ", length(slots3[2]))
println("Heavy machinery: ", length(slots3[3]))
println("Secu: ", length(slots3[4]))

deck4 = filter(x -> x.deck_id == 4, slots)
slots4 = [filter(x->x.cargo_type_id == i, deck4) for i in 1:4]
println("Deck 4")
println("trucks: ", length(slots4[1]))
println("cars: ", length(slots4[2]))
println("Heavy machinery: ", length(slots4[3]))
println("Secu: ", length(slots4[4]))

deck5 = filter(x -> x.deck_id == 5, slots)
slots5 = [filter(x -> x.cargo_type_id == i, deck5) for i in 1:4]
println("Deck 5")
println("trucks: ", length(slots5[1]))
println("cars: ", length(slots5[2]))
println("Heavy machinery: ", length(slots5[3]))
println("Secu: ", length(slots5[4]))
=#

for i in 1:length(Hollandia_test)
    test_instance_hol = Hollandia_test[i]
    test_instance_fin = Finlandia_test[i]
Deterministic_problem_hol = load_data("hollandia",test_instance_hol,"hazardous")
Deterministic_problem_fin = load_data("finlandia",test_instance_fin,"hazardous")

    println("Finlandia test $(i)")
    println("Total weight: ", Deterministic_problem_fin.cargo.total_weight)
    println("number of cars:", length(filter(x -> x.cargo_type_id == 2, Deterministic_problem_fin.cargo)))
    println("number of trucks:", length(filter(x -> x.cargo_type_id == 1, Deterministic_problem_fin.cargo)))
    println("number of Secu-boxes:", length(filter(x -> x.cargo_type_id == 4, Deterministic_problem_fin.cargo)))
    println("number of heavy machinery:", length(filter(x -> x.cargo_type_id == 3, Deterministic_problem_fin.cargo)))
    
    println("Hollandia test $(i)")
    println("Total weight: ", Deterministic_problem_hol.cargo.total_weight)
    println("number of cars:", length(filter(x -> x.cargo_type_id == 2, Deterministic_problem_hol.cargo)))
    println("number of trucks:", length(filter(x -> x.cargo_type_id == 1, Deterministic_problem_hol.cargo)))
    println("number of Secu-boxes:", length(filter(x -> x.cargo_type_id == 4, Deterministic_problem_hol.cargo)))
    println("number of heavy machinery:", length(filter(x -> x.cargo_type_id == 3, Deterministic_problem_hol.cargo))) 
end