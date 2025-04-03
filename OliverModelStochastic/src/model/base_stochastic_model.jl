# Base stochastic model


function create_model_stochastic(problem::StochasticStowageProblem)

	@unpack vessel, slots, cargo, unknown_weights, known_weights, scenarios, probability = problem

	n_slots = length(slots)
	cargo_types = cargo.items[1].cargo_types # types are the same across all scenarios
	n_positions = length(vessel.frame_positions)
	n_cargo = length(cargo.items[1])
	n_deck = length(vessel.decks)
	n_ballast_tanks = length(vessel.ballast_tanks)
    Cargo_1_scenario = cargo.items[1]

	# Rasmus: Fraction of slot within this frame section
	slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)

    # Water density
	ρ = 1.025

	#model = Model(() -> HiGHS.Optimizer())
	#set_time_limit_sec(model, 60 * 0.5)
	#model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
	model = Model(Gurobi.Optimizer)
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4) 
	# model = Model(HiGHS.Optimizer)

	# Weight variables
	@variable(model, weight[1:n_slots,1:scenarios] >= 0)   # Weight at each slot, in each scenario

	# Cargo assignment
	# Rasmus: cs is a binary matrix, 1 if cargo c is assigned to slot s
	# is called a_{c,s} in paper
	@variable(model, cs[1:n_cargo, 1:n_slots], Bin)  # Assignment variables
	# Rasmus: cargo_slack is x_c in paper
	@variable(model, cargo_slack[1:n_cargo], Bin) # 1 if cargo is assigned a slot



	function check_refrigeration_compatibility(cargo::Cargo, slot::Slot)::Bool
		# If cargo needs refrigeration, slot must be refrigerated
		if cargo.refers
			return slot.refrigerated
		end
		# Non-refrigerated cargo can go anywhere
		return true
	end

    # Since refrigeration is constant in all scenarios, this doesn't need to change
	for c in eachindex(cargo.items[1])
		for s in eachindex(slots)
			# If assignment is incompatible, force x[i,j] to be 0
			if !check_refrigeration_compatibility(Cargo_1_scenario[c], slots[s])
				fix(cs[c, s], 0; force = true)
			end
		end
	end

	# This is a bit confusing
	CSC = sum(vessel.ballast_tanks[t].max_vol for t in 1:n_ballast_tanks)
	haz_cargo = [c.hazardous for c in Cargo_1_scenario] .!= 18 # Rasmus: Boolean array, why 18?
	cost = [CSC for c in Cargo_1_scenario]
	cost = cost .+ haz_cargo * CSC # Rasmus: Pretty sure this is the pseudo-revenue for the objective function
	area = [get_length(Cargo_1_scenario[c]) * get_width(Cargo_1_scenario[c]) for c in 1:n_cargo]

	# One cargo per slot at most. Constraint (24)
	@constraint(model, [s = 1:n_slots],
		sum(cs[c, s] for c ∈ 1:n_cargo) <= 1)
	# Constraint (25)
	@constraint(model,[c = 1:n_cargo], 
	sum(cs[c, s] for s ∈ 1:n_slots) <= 1)

	# Slot weight calculation. Constraint (27)
	@constraint(model, [s = 1:n_slots, sc = 1:scenarios],
        weight[s,sc] == sum(cargo.items[sc][c].weight * cs[c, s] for c ∈ 1:n_cargo))
		#weight[s] == sum(cargo[c].weight * cs[c, s] for c ∈ 1:n_cargo))

	# Defining cargo_slack variable. Constraint (23)
	@constraint(model,[c = 1:n_cargo],
		sum(cs[c, s] for s ∈ 1:n_slots) == cargo_slack[c]
	)
	# Rasmus: Function creates a list of slot ids which are not equal to the cargo type id input "t"
	invalid_slots(t::Int) = [slot.id for slot in filter(x -> x.cargo_type_id != t, slots)]
	# Rasmus: Sets all slots which are not compatible with the cargo type to 0
	@constraint(model,
		[t in cargo_types,
			i in [cargo.id for cargo in filter(x -> x.cargo_type_id == t, Cargo_1_scenario)]],
		sum(cs[i, j] for j in invalid_slots(t)) == 0)
	# Rasmus: Overlapping slots 
	overlapping_indicies = [Tuple(ix) for ix in findall(slots.overlap_matrix)]
	# Rasmus: Can only use one of two slots if they overlap. Constraint (26)
	@constraint(model, [(x, y) in overlapping_indicies], sum(cs[:, x]) + sum(cs[:, y]) <= 1)

	# Rasmus: expression defines something (pos_weight_cargo) without it being a variable or constraint
	# Rasmus: pos_weight_cargo is used as input later
	# Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions, sc = 1:scenarios],
		# Cargo weights using precalculated proportions
		sum(weight[s,sc] * slots_to_frame[s, p] for s ∈ 1:n_slots)
	)
	# Deck weight limit. Constraint (29)
	@constraint(model,[d = 1:n_deck, sc = 1:scenarios],
		sum(weight[s,sc] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit)
	# Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
	@expression(model, lcg_cargo[sc = 1:scenarios], sum(weight[s,sc] * slots[s].lcg for s ∈ 1:n_slots))
	@expression(model, tcg_cargo[sc = 1:scenarios], sum(weight[s,sc] * slots[s].tcg for s ∈ 1:n_slots))
	# Rasmus: Don't understand why is not using weight[s]
	@expression(model, vcg_cargo[sc = 1:scenarios], sum(cs[c, s] * cargo.items[sc].items[c].weight * (slots[s].vcg + get_height(cargo.items[1][c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

	# Rasmus: y_t is defined in this part of the model
	add_stability_stochastic!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo, scenarios)

	ballast_volume = model[:ballast_volume]

	@objective(model, Min,
		sum(ballast_volume[t,sc]*probability[sc] for t ∈ 1:n_ballast_tanks, sc in 1:scenarios)
		-
		sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
	)

	return model
end

# if length([cargo.haz_class for cargo in cargos if cargo.haz_class != 18]) > 0
#     @info "Adding hazardous cargo constraints"
#     add_hazardous_cargo_constraints!(model, vessel, cargos, m)
# end



# optimize!(model)

# solution = extract_solution(problem, model)

# plot_cargo_solution(solution::Solution)
