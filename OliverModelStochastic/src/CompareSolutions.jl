# Compare solutions

function compare_solutions_print(stochastic_solution::SolutionStochastic, deterministic_solution::Solution,sc::Int64 = nothing)#problem::StowageProblem)
    println("Total number of cargo: ", stochastic_solution.n_cargo_total)
    println("--------------------------------------")
    println("Scenario(s): ", sc == nothing ? "All" : sc)
    println("Status of stochastic model: ", stochastic_solution.status)
    println("Stochastic objective value: ", stochastic_solution.objective)
    println("Stochastic model time: ", stochastic_solution.time)
    println("Stochastic model cargo loaded: ", stochastic_solution.n_cargo_loaded)
    println("Number of scenarios: ", length(stochastic_solution.forces))
    println("--------------------------------------")
    if !isnothing(sc)
        println("Objective value for given scenario: ", stochastic_solution.objective_scenario[sc])
        println("Printing forces for scenario: ", sc)
        println("Total cargo weight: ", stochastic_solution.forces[sc].cargo_weight)
        println("Total weight: ", stochastic_solution.forces[sc].total_weight)
        println("Ballast weight: ", stochastic_solution.forces[sc].ballast_weight)
        #println("Shear force: ", stochastic_solution.forces[sc].shear_force)
        #println("Bending moment: ", stochastic_solution.forces[sc].bending_moment)
        #println("Ballast volume: ", stochastic_solution.forces[sc].ballast_volume)
        println("LCG: ", stochastic_solution.forces[sc].lcg)
        println("TCG: ", stochastic_solution.forces[sc].tcg)
        println("VCG: ", stochastic_solution.forces[sc].vcg)
    end
    println("--------------------------------------")
    println("Status of deterministic model: ", deterministic_solution.status)
    println("Deterministic objective value: ", deterministic_solution.objective)
    println("Deterministic model time: ", deterministic_solution.time)
    println("Deterministic model cargo loaded: ", deterministic_solution.n_cargo_loaded)
    println("Printing forces for deterministic scenario:")
    println("Total cargo weight: ", deterministic_solution.cargo_weight)
    println("Total weight: ", deterministic_solution.total_weight)
    println("Ballast weight: ", deterministic_solution.ballast_weight)
    #println("Shear force: ", deterministic_solution.shear_force)
    #println("Bending moment: ", deterministic_solution.bending_moment)
    #println("Ballast volume: ", deterministic_solution.ballast_volume)
    println("LCG: ", deterministic_solution.lcg)
    println("TCG: ", deterministic_solution.tcg)
    println("VCG: ", deterministic_solution.vcg)
    println("--------------------------------------")
    println("Model statistics: ")
    println("Stochastic model")
    println("Number of variables: ", stochastic_solution.n_variables)
    println("Number of constraints: ", stochastic_solution.n_constraints)
    println("Size: $(round(stochastic_solution.model_size/1024/1024, digits=2)) MB")
    println("Deterministic model")
    println("Number of variables: ", deterministic_solution.n_variables)
    println("Number of constraints: ", deterministic_solution.n_constraints)
    println("Size: $(round(deterministic_solution.model_size/1024/1024, digits=2)) MB")
    println("--------------------------------------")
    diff_ballast_weight = Float64[]
    diff_cargo_weight = Float64[]
    for i in 1:length(stochastic_solution.forces)
        push!(diff_ballast_weight, stochastic_solution.forces[i].ballast_weight - deterministic_solution.ballast_weight)
        push!(diff_cargo_weight, stochastic_solution.forces[i].cargo_weight - deterministic_solution.cargo_weight)
        println("Scenario: ",i)
        println("Difference in ballast_weight: ", stochastic_solution.forces[i].ballast_weight - deterministic_solution.ballast_weight)
    end
    diff, diff_index, same_index = compare_placements(stochastic_solution, deterministic_solution)
    println("Number of different placements: ", diff)
    println("Indices of different placements: ", diff_index)
    println("Indices of same placements: ", same_index)
    println("--------------------------------------")

    return diff_ballast_weight, diff_cargo_weight
end
# Compare the placements of cargo for stochastic and deterministic model
function compare_placements(stochastic_solution::SolutionStochastic, deterministic_solution::Solution,)
    diff = 0
    diff_index = []
    same_index = []
    for i in 1:minimum([stochastic_solution.n_cargo_loaded, deterministic_solution.n_cargo_loaded])
        if !same_placement(stochastic_solution.cargo[i], deterministic_solution.cargo[i])
            diff += 1
            push!(diff_index, stochastic_solution.cargo[i].id)
        else
            push!(same_index, stochastic_solution.cargo[i].id)
        end
    end
    return diff, diff_index, same_index
end
function same_placement(cargo1, cargo2)
    return ((cargo1.lcg == cargo2.lcg) && (cargo1.tcg == cargo2.tcg) && (cargo1.vcg == cargo2.vcg))
end
