# Given placement and realization, optimize ballast water 

# cs is the placement for the cargo, decided in the first stage model.
# problem is the deterministic version of the problem.
function second_stage_model(cs, problem::StowageProblem)

    @unpack vessel, slots, cargo = problem

	n_slots = length(slots)
	cargo_types = cargo.cargo_types
	n_positions = length(vessel.frame_positions)
	n_cargo = length(cargo)
	n_deck = length(vessel.decks)
	n_ballast_tanks = length(vessel.ballast_tanks)
    
    slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)
    ρ = 1.025
	model = Model(Gurobi.Optimizer)
    @variable(model, weight[1:n_slots] >= 0)   # Weight at each slot
    
    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
	haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
	cost = [CSC for c in cargo]
	cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
	area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]

    # Slot weight calculation. Constraint (27)
	@constraint(model, [s = 1:n_slots],
    weight[s] == sum(cargo[c].weight * cs[c, s] for c ∈ 1:n_cargo))

    # Expression for cargo_slack
    @expression(model,cargo_slack[c = 1:n_cargo],
		sum(cs[c, s] for s ∈ 1:n_slots))

    # Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions],
    # Cargo weights using precalculated proportions
    sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots)
)
    # # # Deck weight limit. Constraint (29)
	@constraint(model,[d = 1:n_deck],
		sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit
	)
    


    # Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
	@expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
	@expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
	# Rasmus: Don't understand why is not using weight[s]
	@expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

    add_stability!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)

	ballast_volume = model[:ballast_volume]

	@objective(model, Min,
		sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
		-
		sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
	)
	return model
end