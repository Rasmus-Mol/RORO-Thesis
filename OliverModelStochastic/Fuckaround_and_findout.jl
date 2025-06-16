include("packages_and_files.jl")

# Do stuff 
test_instance_hol = Hollandia_test[5]
test_instance_fin = Finlandia_test[1]

Deterministic_problem_hol = load_data("hollandia",test_instance_hol,"hazardous")
#Deterministic_problem_fin = load_data("finlandia","no_cars_medium_100_haz_eq_0.1","hazardous")
Deterministic_problem_fin = load_data("finlandia",test_instance_fin,"hazardous")

ids = [c.id for c in filter(x-> x.cargo_type_id != 2, Deterministic_problem_fin.cargo)]
ids_cars = [c.id for c in filter(x-> x.cargo_type_id == 2, Deterministic_problem_fin.cargo)]

temp = create_stochastic_problem_cars_known(Deterministic_problem_fin, 10)

# 
model = Model(Gurobi.Optimizer)
# number should match number of cores used at HPC
set_optimizer_attribute(model, "Threads", 4)
@variable(model, x_1 >= 0)
@objective(model, Min, x_1)
# Solve the model
optimize!(model)
a = primal_status(model)
ts = termination_status(model)
if ts in (MOI.OPTIMAL, MOI.FEASIBLE_POINT)
    println(2)
end
if a == MOI.FEASIBLE_POINT
    println("Feasible point found")
end
import MathOptInterface as MOI
#
Deterministic_problem_fin.vessel.decks
slots_deck1 = filter(x -> x.deck_id == 1, Deterministic_problem_fin.slots)
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
pro_fin_1 = create_stochastic_problem(Deterministic_problem_fin, 100, length(Deterministic_problem_fin.cargo), [], Bootstrap_bookedweight_quantile)

pro_fin_2 = create_stochastic_problem_scenarioreduction(Deterministic_problem_fin, 10,
length(Deterministic_problem_fin.cargo), 100, [], Bootstrap_bookedweight_quantile, scenario_reduction_heuristic, 5*60)

pro_fin_3 = create_stochastic_problem_scenarioreduction(Deterministic_problem_fin, 10,
length(Deterministic_problem_fin.cargo), 50, [], Bootstrap_bookedweight_quantile, scenario_reduction_naive, 60)


initial_scenarios = 100
problem_det_noise = add_white_noise_to_test_instance(Deterministic_problem_fin)
pro = create_stochastic_problem(problem_det_noise, initial_scenarios, length(Deterministic_problem_fin.cargo), [],Bootstrap_bookedweight_quantile) 
pro_reduced = scenario_reduced(pro, 10, scenario_reduction_clustering,60)
pro.probability
println(pro_reduced.probability)

c,p = scenario_reduction_clustering(pro.cargo, pro.probability, 10, 60)
p
CargoC = pro.cargo
sc_old = length(CargoC.items)
    cost = zeros(sc_old, sc_old)
    for i in 1:sc_old
        for j in 1:sc_old
            cost[i,j] = sum(abs.(getfield.(CargoC.items[i],:weight) .- getfield.(CargoC.items[j],:weight)))
        end
    end
    heatmap(cost)
    # Perform hierarchical clustering
    result = hclust(cost, linkage = :ward)  # or :single, :complete, :ward, etc.
    # Cut the dendrogram into sc clusters
    labels = cutree(result, k = 10)
    if sum([count(x->x == i, labels) for i in 1:10]) != sc_old
        throw("Clustering did not produce the expected number of scenarios")
    end
    println(labels)
    av = [count(x->x == i, labels) for i in 1:10]
    wa = [count(x->x == i, labels) for i in 1:10]
    com = [count(x->x == i, labels) for i in 1:10]
    println(av)
    println(wa)
    println(com)


pro_biased_noise = add_biased_noise_to_test_instance(Deterministic_problem_fin)
println(mean([c.weight for c in filter(x->x.cargo_type_id == 1, pro_biased_noise.cargo.items)]))
println(mean([c.weight for c in filter(x->x.cargo_type_id == 1, Deterministic_problem_fin.cargo.items)]))
println( mean([mean([c.weight for c in filter(x->x.cargo_type_id == 1, pro_fin_1.cargo.items[i].items)]) 
for i in pro_fin_1.scenarios]))
println(mean([c.weight for c in filter(x->x.cargo_type_id == 4, pro_biased_noise.cargo.items)]))
println(mean([c.weight for c in filter(x->x.cargo_type_id == 4, Deterministic_problem_fin.cargo.items)]))
println( mean([mean([c.weight for c in filter(x->x.cargo_type_id == 4, pro_fin_1.cargo.items[i].items)]) 
for i in pro_fin_1.scenarios]))



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

CargoC = pro_fin.cargo
probability = pro_fin.probability
lab, CC, prob = scenario_reduction_clustering(CargoC, probability, 10, 5)
sum(prob)
println([count(x->x == i, lab) for i in 1:10])
println(CC.items.total_weight)
CargoC.items.total_weight
idx = findall(t -> t == 2, lab)
cargo_scenarios = Vector{CargoCollection}()
prob_new = zeros(length(idx))

test = scenario_reduced(pro_fin, 10, 
scenario_reduction_clustering,60)
test.probability

println(lab)


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

