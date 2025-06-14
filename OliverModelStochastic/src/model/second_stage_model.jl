# Given placement and realization, optimize ballast water 

# cs is the placement for the cargo, decided in the first stage model.
# problem is the deterministic version of the problem.
function second_stage_model(cs1, problem::StowageProblem)

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
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)
    @variable(model, weight[1:n_slots] >= 0)   # Weight at each slot
    @variable(model, cs[1:n_cargo, 1:n_slots], Bin)  # Assignment variables
	@variable(model, cargo_slack[1:n_cargo], Bin) # 1 if cargo is assigned a slot
	# Ensure placement is the same
	@constraint(model, [c = 1:n_cargo, s = 1:n_slots],
		cs[c, s] == cs1[c, s])


    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
	haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
	cost = [CSC for c in cargo]
	cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
	area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]

	# One cargo per slot at most. Constraint (24)
	@constraint(model, [s = 1:n_slots],
		sum(cs[c, s] for c ∈ 1:n_cargo) <= 1)
	# One slot per cargo at most. Constraint not in paper. Constraint (25)
	@constraint(
		model,
		[c = 1:n_cargo], sum(cs[c, s] for s ∈ 1:n_slots) <= 1)


    # Slot weight calculation. Constraint (27)
	@constraint(model, [s = 1:n_slots],
    weight[s] == sum(cargo[c].weight * cs[c, s] for c ∈ 1:n_cargo))

    # Expression for cargo_slack
    #@expression(model,cargo_slack[c = 1:n_cargo],
	#	sum(cs[c, s] for s ∈ 1:n_slots))
	@constraint(model,[c = 1:n_cargo],
		sum(cs[c, s] for s ∈ 1:n_slots) == cargo_slack[c]
	)
	# Rasmus: Function creates a list of slot ids which are not equal to the cargo type id input "t"
	invalid_slots(t::Int) = [slot.id for slot in filter(x -> x.cargo_type_id != t, slots)]
	# Rasmus: Sets all slots which are not compatible with the cargo type to 0
	@constraint(model,
		[t in cargo_types,
			i in [cargo.id for cargo in filter(x -> x.cargo_type_id == t, cargo)]],
		sum(cs[i, j] for j in invalid_slots(t)) == 0)
	# Rasmus: Overlapping slots 
	overlapping_indicies = [Tuple(ix) for ix in findall(slots.overlap_matrix)]
	# Rasmus: Can only use one of two slots if they overlap. Constraint (26)
	@constraint(model, [(x, y) in overlapping_indicies], sum(cs[:, x]) + sum(cs[:, y]) <= 1)


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

# Same model, just with slack variable to see why the model is infeasible
function second_stage_model_slack(cs1, problem::StowageProblem)
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
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)
    @variable(model, weight[1:n_slots] >= 0)   # Weight at each slot
    @variable(model, cs[1:n_cargo, 1:n_slots], Bin)  # Assignment variables
	@variable(model, cargo_slack[1:n_cargo], Bin) # 1 if cargo is assigned a slot
	# Ensure placement is the same
	@constraint(model, [c = 1:n_cargo, s = 1:n_slots],
		cs[c, s] == cs1[c, s])

    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
	haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
	cost = [CSC for c in cargo]
	cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
	area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]
	
	# Penalty for violating constraints
	M = 100000 # Should be determined more precisely at some point
	#M = sum(cost) +1 

    # Slot weight calculation. Constraint (27)
	@constraint(model, [s = 1:n_slots],
    weight[s] == sum(cargo[c].weight * cs[c, s] for c ∈ 1:n_cargo))

    # Expression for cargo_slack
    #@expression(model,cargo_slack[c = 1:n_cargo],
	#	sum(cs[c, s] for s ∈ 1:n_slots))
	@constraint(model,[c = 1:n_cargo],
		sum(cs[c, s] for s ∈ 1:n_slots) == cargo_slack[c]
	)

	# Rasmus: Function creates a list of slot ids which are not equal to the cargo type id input "t"
	invalid_slots(t::Int) = [slot.id for slot in filter(x -> x.cargo_type_id != t, slots)]
	# Rasmus: Sets all slots which are not compatible with the cargo type to 0
	@constraint(model,
		[t in cargo_types,
			i in [cargo.id for cargo in filter(x -> x.cargo_type_id == t, cargo)]],
		sum(cs[i, j] for j in invalid_slots(t)) == 0)
	# Rasmus: Overlapping slots 
	overlapping_indicies = [Tuple(ix) for ix in findall(slots.overlap_matrix)]
	# Rasmus: Can only use one of two slots if they overlap. Constraint (26)
	@constraint(model, [(x, y) in overlapping_indicies], sum(cs[:, x]) + sum(cs[:, y]) <= 1)

    # Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions],
    # Cargo weights using precalculated proportions
    sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots))
    # # # Deck weight limit. Constraint (29)
	@variable(model, slack_deck[1:n_deck] >= 0)
	@constraint(model,[d = 1:n_deck],
		sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit + slack_deck[d]
	)
    
    # Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
	@expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
	@expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
	# Rasmus: Don't understand why is not using weight[s]
	@expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

	# Slack variables:
	# Center of gravity
	@variable(model, slack_Vmax >=0)
	@variable(model, slack_Vmin >=0)
	@variable(model, slack_Tmin>=0)
	@variable(model, slack_Tmax>=0)
	@variable(model, slack_Lmin>=0)
	@variable(model, slack_Lmax>=0)
	# Stress and Bending
	@variable(model, slack_shear1[1:n_positions]>=0)
	@variable(model, slack_shear2[1:n_positions]>=0)
	stress_limits = vessel.stress_limits
	@variable(model, slack_shearMin[1:length(stress_limits)]>=0)
	@variable(model, slack_shearMax[1:length(stress_limits)]>=0)
	@variable(model, slack_bendingMax[1:length(stress_limits)]>=0)
	# Ballast tanks
	@variable(model, slack_ballast_tanks[1:n_ballast_tanks] >= 0)

    add_stability_slack!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo,
	slack_Vmax,slack_Vmin,slack_Tmin,slack_Tmax,slack_Lmin,slack_Lmax,
	slack_shear1,slack_shear2,slack_shearMin,slack_shearMax,slack_bendingMax,slack_ballast_tanks)

	ballast_volume = model[:ballast_volume]

	@objective(model, Min,
		sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
		-
		sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
		+ M * (slack_Lmax+slack_Lmin+slack_Tmax+slack_Tmin+slack_Vmax+slack_Vmin
		+ sum(slack_shear1) + sum(slack_shear2)+sum(slack_shearMax) + sum(slack_shearMin) + sum(slack_bendingMax)
		+ sum(slack_deck)+sum(slack_ballast_tanks))
		)
	return model
end

# Slack on everything, initial model had slack on deck weight only.
function second_stage_model_deckslacked_slack(cs1, problem::StowageProblem,slack_fraction::Vector{Float64})
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
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)
    @variable(model, weight[1:n_slots] >= 0)   # Weight at each slot
    @variable(model, cs[1:n_cargo, 1:n_slots], Bin)  # Assignment variables
	@variable(model, cargo_slack[1:n_cargo], Bin) # 1 if cargo is assigned a slot
	# Ensure placement is the same
	@constraint(model, [c = 1:n_cargo, s = 1:n_slots],
		cs[c, s] == cs1[c, s])

    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
	haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
	cost = [CSC for c in cargo]
	cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
	area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]
	
	# Penalty for violating constraints
	M = 100000 # Should be determined more precisely at some point
	#M = sum(cost) +1 

    # Slot weight calculation. Constraint (27)
	@constraint(model, [s = 1:n_slots],
    weight[s] == sum(cargo[c].weight * cs[c, s] for c ∈ 1:n_cargo))

    # Expression for cargo_slack
    #@expression(model,cargo_slack[c = 1:n_cargo],
	#	sum(cs[c, s] for s ∈ 1:n_slots))
	@constraint(model,[c = 1:n_cargo],
		sum(cs[c, s] for s ∈ 1:n_slots) == cargo_slack[c]
	)

	# Rasmus: Function creates a list of slot ids which are not equal to the cargo type id input "t"
	invalid_slots(t::Int) = [slot.id for slot in filter(x -> x.cargo_type_id != t, slots)]
	# Rasmus: Sets all slots which are not compatible with the cargo type to 0
	@constraint(model,
		[t in cargo_types,
			i in [cargo.id for cargo in filter(x -> x.cargo_type_id == t, cargo)]],
		sum(cs[i, j] for j in invalid_slots(t)) == 0)
	# Rasmus: Overlapping slots 
	overlapping_indicies = [Tuple(ix) for ix in findall(slots.overlap_matrix)]
	# Rasmus: Can only use one of two slots if they overlap. Constraint (26)
	@constraint(model, [(x, y) in overlapping_indicies], sum(cs[:, x]) + sum(cs[:, y]) <= 1)

    # Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions],
    # Cargo weights using precalculated proportions
    sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots))
    # # # Deck weight limit. Constraint (29)
	@variable(model, slack_deck[1:n_deck] >= 0)
	@constraint(model,[d = 1:n_deck],
		sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit * slack_fraction[d] + slack_deck[d]
	)
    
    # Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
	@expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
	@expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
	# Rasmus: Don't understand why is not using weight[s]
	@expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

	# Slack variables:
	# Center of gravity
	@variable(model, slack_Vmax >=0)
	@variable(model, slack_Vmin >=0)
	@variable(model, slack_Tmin>=0)
	@variable(model, slack_Tmax>=0)
	@variable(model, slack_Lmin>=0)
	@variable(model, slack_Lmax>=0)
	# Stress and Bending
	@variable(model, slack_shear1[1:n_positions]>=0)
	@variable(model, slack_shear2[1:n_positions]>=0)
	stress_limits = vessel.stress_limits
	@variable(model, slack_shearMin[1:length(stress_limits)]>=0)
	@variable(model, slack_shearMax[1:length(stress_limits)]>=0)
	@variable(model, slack_bendingMax[1:length(stress_limits)]>=0)
	# Ballast tanks
	@variable(model, slack_ballast_tanks[1:n_ballast_tanks] >= 0)

    add_stability_slack!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo,
	slack_Vmax,slack_Vmin,slack_Tmin,slack_Tmax,slack_Lmin,slack_Lmax,
	slack_shear1,slack_shear2,slack_shearMin,slack_shearMax,slack_bendingMax,slack_ballast_tanks)

	ballast_volume = model[:ballast_volume]

	@objective(model, Min,
		sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
		-
		sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
		+ M * (slack_Lmax+slack_Lmin+slack_Tmax+slack_Tmin+slack_Vmax+slack_Vmin
		+ sum(slack_shear1) + sum(slack_shear2)+sum(slack_shearMax) + sum(slack_shearMin) + sum(slack_bendingMax)
		+ sum(slack_deck)+sum(slack_ballast_tanks))
		)
	return model
end
# Second stage model with only slack on decks
function second_stage_model_deckslacked(cs1, problem::StowageProblem, slack_fraction::Vector{Float64})

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
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)
    @variable(model, weight[1:n_slots] >= 0)   # Weight at each slot
    @variable(model, cs[1:n_cargo, 1:n_slots], Bin)  # Assignment variables
	@variable(model, cargo_slack[1:n_cargo], Bin) # 1 if cargo is assigned a slot
	# Ensure placement is the same
	@constraint(model, [c = 1:n_cargo, s = 1:n_slots],
		cs[c, s] == cs1[c, s])


    CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
	haz_cargo = [c.hazardous for c in cargo] .!= 18 # Rasmus: Boolean array, why 18?
	cost = [CSC for c in cargo]
	cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
	area = [get_length(cargo[c]) * get_width(cargo[c]) for c in 1:n_cargo]

	# One cargo per slot at most. Constraint (24)
	@constraint(model, [s = 1:n_slots],
		sum(cs[c, s] for c ∈ 1:n_cargo) <= 1)
	# One slot per cargo at most. Constraint not in paper. Constraint (25)
	@constraint(
		model,
		[c = 1:n_cargo], sum(cs[c, s] for s ∈ 1:n_slots) <= 1)


    # Slot weight calculation. Constraint (27)
	@constraint(model, [s = 1:n_slots],
    weight[s] == sum(cargo[c].weight * cs[c, s] for c ∈ 1:n_cargo))

    # Expression for cargo_slack
    #@expression(model,cargo_slack[c = 1:n_cargo],
	#	sum(cs[c, s] for s ∈ 1:n_slots))
	@constraint(model,[c = 1:n_cargo],
		sum(cs[c, s] for s ∈ 1:n_slots) == cargo_slack[c]
	)
	# Rasmus: Function creates a list of slot ids which are not equal to the cargo type id input "t"
	invalid_slots(t::Int) = [slot.id for slot in filter(x -> x.cargo_type_id != t, slots)]
	# Rasmus: Sets all slots which are not compatible with the cargo type to 0
	@constraint(model,
		[t in cargo_types,
			i in [cargo.id for cargo in filter(x -> x.cargo_type_id == t, cargo)]],
		sum(cs[i, j] for j in invalid_slots(t)) == 0)
	# Rasmus: Overlapping slots 
	overlapping_indicies = [Tuple(ix) for ix in findall(slots.overlap_matrix)]
	# Rasmus: Can only use one of two slots if they overlap. Constraint (26)
	@constraint(model, [(x, y) in overlapping_indicies], sum(cs[:, x]) + sum(cs[:, y]) <= 1)


    # Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions],
    # Cargo weights using precalculated proportions
    sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots)
	)
    # # # Deck weight limit. Constraint (29)
	@constraint(model,[d = 1:n_deck],
		sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit * slack_fraction[d]
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
