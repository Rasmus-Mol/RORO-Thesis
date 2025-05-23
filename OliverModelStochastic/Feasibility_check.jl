# Checks if solution can be feasible if we change weights of cargo
function feasibility_check(sol::Solution, problem::StowageProblem, CargoC::CargoCollection)
    cs = sol.cs
    #tanks = Solution.ballast_volume
    slots = problem.slots
    vessel = problem.vessel
    name = problem.name
    sz = size(cs)
    n_slots = sz[2]
	cargo_types = CargoC.cargo_types
	n_positions = length(vessel.frame_positions)
	n_cargo = sz[1]
    if n_cargo != length(CargoC)
        error("Number of cargo in cs and CargoC do not match")
    end
	n_deck = length(vessel.decks)
	n_ballast_tanks = length(vessel.ballast_tanks)
    # Feasibility:
    feasible = true
    # Water density.
    ρ = 1.025
    # Weight in slots
    weight = zeros(n_slots)
    weight[1:n_slots] = [sum(CargoC[c].weight * cs[c, s] for c ∈ 1:n_cargo) for s in 1:n_slots]
    # Deck weight
    for d in 1:n_deck
        deck_weight = sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)])
        if deck_weight > vessel.decks[d].weight_limit
            str = "Deck weight limit exceeded on deck $(d)"
            println(str)
            feasible = false
            return feasible, str, 0
        end
    end
    # Optimize stability with locked cargo placements and new cargo weights.
    pro = StowageProblem(
        # Core components
    vessel = vessel,
    slots = slots,
    cargo = CargoC,  
    name = name,
    timestamp = now()
    )
    # Create and optimize model
    model = second_stage_model_v2(cs, pro)
    set_time_limit_sec(model, 60 * 5) # 5 minutes
    set_silent(model)
    optimize!(model)
    # Check if the model is feasible
    if termination_status(model) == MOI.INFEASIBLE
        str = "Model is infeasible - cannot be stable with this weight."
        println(str)
        feasible = false
        return feasible, str, 0
    end
    if termination_status(model) == MOI.TIMEOUT
        str = "Stability Model timed out"
        println(str)
        feasible = false
        return feasible, str, model
    end
    if termination_status(model) == MOI.OPTIMAL
        str = "Stability Model is feasible"
        println(str)
        feasible = true
        return feasible, str, model
    end
end

