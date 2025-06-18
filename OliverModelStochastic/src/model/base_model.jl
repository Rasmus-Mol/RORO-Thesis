
function create_model(problem::StowageProblem)

	@unpack vessel, slots, cargo = problem

	n_slots = length(slots)
	cargo_types = cargo.cargo_types
	n_positions = length(vessel.frame_positions)
	n_cargo = length(cargo)
	n_deck = length(vessel.decks)
	n_ballast_tanks = length(vessel.ballast_tanks)

	# Rasmus: Fraction of slot within this frame section
	slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)

	# draft_index_min, 
	# draft_index_max,
	# n_positions,

	# Unpack all needed values from common_params
	# @unpack n_ballast_tanks,
	#         ballast_tank_weight_max,
	#         ballast_tank_proportions,
	#         ship_position_weight,
	#         ship_buoyancy_weight,
	#         ship_buoyancy_position,
	#         force_index,
	#         shear_min,
	#         shear_max,
	#         bending_max = problem.vessel

	# Unpack needed values from solver struct m
	# @unpack parameter, ship = m
	# @unpack CSC, D, draft_index_full, TCGmax, TCGmin, vcgnum, lcgnum, B = parameter
	# @unpack regulartank, hydro, deck = ship
	# @unpack lcg, tcg = regulartank
	# @unpack displacement, lcb, gm_min, metacenter, kg_min = hydro

	ρ = 1.025

	#model = Model(() -> HiGHS.Optimizer())
	#set_time_limit_sec(model, 60 * 0.5)
	# model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
	model = Model(Gurobi.Optimizer)
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)
	# model = Model(HiGHS.Optimizer)

	# Weight variables
	@variable(model, weight[1:n_slots] >= 0)   # Weight at each slot

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

	for c in eachindex(cargo)
		for s in eachindex(slots)
			# If assignment is incompatible, force x[i,j] to be 0
			if !check_refrigeration_compatibility(cargo[c], slots[s])
				fix(cs[c, s], 0; force = true)
			end
		end
	end

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

	# Defining cargo_slack variable. Constraint (23)
	@constraint(
		model,
		[c = 1:n_cargo],
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

	# Rasmus: expression defines something (pos_weight_cargo) without it being a variable or constraint
	# Rasmus: pos_weight_cargo is used as input later
	# Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions],
		# Cargo weights using precalculated proportions
		sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots)
	)

	# # # Deck weight limit. Constraint (29)
	@constraint(
		model,
		[d = 1:n_deck],
		sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit
	)

	# Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
	@expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
	@expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
	# Rasmus: Don't understand why is not using weight[s]
	@expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

	# Rasmus: y_t is defined in this part of the model
	add_stability!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)

	ballast_volume = model[:ballast_volume]

	@objective(model, Min,
		sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
		-
		sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
	)

	return model
end

# Slacking weight limit on decks
function create_model_slack_deck(problem::StowageProblem,slack_fraction::Vector{Float64})

	@unpack vessel, slots, cargo = problem

	n_slots = length(slots)
	cargo_types = cargo.cargo_types
	n_positions = length(vessel.frame_positions)
	n_cargo = length(cargo)
	n_deck = length(vessel.decks)
	n_ballast_tanks = length(vessel.ballast_tanks)

	# Rasmus: Fraction of slot within this frame section
	slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)

	# draft_index_min, 
	# draft_index_max,
	# n_positions,

	# Unpack all needed values from common_params
	# @unpack n_ballast_tanks,
	#         ballast_tank_weight_max,
	#         ballast_tank_proportions,
	#         ship_position_weight,
	#         ship_buoyancy_weight,
	#         ship_buoyancy_position,
	#         force_index,
	#         shear_min,
	#         shear_max,
	#         bending_max = problem.vessel

	# Unpack needed values from solver struct m
	# @unpack parameter, ship = m
	# @unpack CSC, D, draft_index_full, TCGmax, TCGmin, vcgnum, lcgnum, B = parameter
	# @unpack regulartank, hydro, deck = ship
	# @unpack lcg, tcg = regulartank
	# @unpack displacement, lcb, gm_min, metacenter, kg_min = hydro

	ρ = 1.025

	#model = Model(() -> HiGHS.Optimizer())
	#set_time_limit_sec(model, 60 * 0.5)
	# model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
	model = Model(Gurobi.Optimizer)
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)
	# model = Model(HiGHS.Optimizer)

	# Weight variables
	@variable(model, weight[1:n_slots] >= 0)   # Weight at each slot

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

	for c in eachindex(cargo)
		for s in eachindex(slots)
			# If assignment is incompatible, force x[i,j] to be 0
			if !check_refrigeration_compatibility(cargo[c], slots[s])
				fix(cs[c, s], 0; force = true)
			end
		end
	end

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

	# Defining cargo_slack variable. Constraint (23)
	@constraint(
		model,
		[c = 1:n_cargo],
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

	# Rasmus: expression defines something (pos_weight_cargo) without it being a variable or constraint
	# Rasmus: pos_weight_cargo is used as input later
	# Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions],
		# Cargo weights using precalculated proportions
		sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots)
	)

	# # # Deck weight limit. Constraint (29)
	@constraint(
		model,
		[d = 1:n_deck],
		sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit * slack_fraction[d]
	)

	# Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
	@expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
	@expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
	# Rasmus: Don't understand why is not using weight[s]
	@expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

	# Rasmus: y_t is defined in this part of the model
	add_stability!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)

	ballast_volume = model[:ballast_volume]

	@objective(model, Min,
		sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
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



# Panic tests
# Olivers Slope
function calculate_vcg_slopes_temp(vessel::Vessel)
	n_tanks = length(vessel.ballast_tanks)
	slopes = zeros(n_tanks)

	for tank_idx in 1:n_tanks
		# Get tank properties
		tank = vessel.ballast_tanks[tank_idx]
		max_vcg = tank.max_vcg
		max_vcg = tank.max_vcg
		min_vcg = tank.min_vcg
		max_vol = tank.max_vol
		#slopes[tank_idx] = round((max_vcg - min_vcg) / max_vol,digits = 2)
		slopes[tank_idx] = max_vcg
	end
	return slopes
end
function create_model_temp(problem::StowageProblem)

	@unpack vessel, slots, cargo = problem

	n_slots = length(slots)
	cargo_types = cargo.cargo_types
	n_positions = length(vessel.frame_positions)
	n_cargo = length(cargo)
	n_deck = length(vessel.decks)
	n_ballast_tanks = length(vessel.ballast_tanks)

	# Rasmus: Fraction of slot within this frame section
	slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)

	# draft_index_min, 
	# draft_index_max,
	# n_positions,

	# Unpack all needed values from common_params
	# @unpack n_ballast_tanks,
	#         ballast_tank_weight_max,
	#         ballast_tank_proportions,
	#         ship_position_weight,
	#         ship_buoyancy_weight,
	#         ship_buoyancy_position,
	#         force_index,
	#         shear_min,
	#         shear_max,
	#         bending_max = problem.vessel

	# Unpack needed values from solver struct m
	# @unpack parameter, ship = m
	# @unpack CSC, D, draft_index_full, TCGmax, TCGmin, vcgnum, lcgnum, B = parameter
	# @unpack regulartank, hydro, deck = ship
	# @unpack lcg, tcg = regulartank
	# @unpack displacement, lcb, gm_min, metacenter, kg_min = hydro

	ρ = 1.025

	#model = Model(() -> HiGHS.Optimizer())
	#set_time_limit_sec(model, 60 * 0.5)
	# model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
	model = Model(Gurobi.Optimizer)
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)
	# model = Model(HiGHS.Optimizer)

	# Weight variables
	@variable(model, weight[1:n_slots] >= 0)   # Weight at each slot

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

	for c in eachindex(cargo)
		for s in eachindex(slots)
			# If assignment is incompatible, force x[i,j] to be 0
			if !check_refrigeration_compatibility(cargo[c], slots[s])
				fix(cs[c, s], 0; force = true)
			end
		end
	end

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

	# Defining cargo_slack variable. Constraint (23)
	@constraint(
		model,
		[c = 1:n_cargo],
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

	# Rasmus: expression defines something (pos_weight_cargo) without it being a variable or constraint
	# Rasmus: pos_weight_cargo is used as input later
	# Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions],
		# Cargo weights using precalculated proportions
		sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots)
	)

	# # # Deck weight limit. Constraint (29)
	@constraint(
		model,
		[d = 1:n_deck],
		sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit
	)

	# Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
	@expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
	@expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
	# Rasmus: Don't understand why is not using weight[s]
	@expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

	# Rasmus: y_t is defined in this part of the model
	add_stability_temp!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)

	ballast_volume = model[:ballast_volume]

	@objective(model, Min,
		sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
		-
		sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
	)

	return model
end
# Rasmus: Add stability constraints to the model
function add_stability_temp!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)
	vcg_slope = calculate_vcg_slopes_temp(vessel)

	n_positions = length(vessel.frame_positions)
	n_ballast_tanks = length(vessel.ballast_tanks)

	# Constants
	ρ = 1.025  # Water density
	TCGmax = 0.001  # TODO: Get from vessel
	TCGmin = -0.001  # TODO: Get from vessel

	# Extract weight dependencies
	draft_index_min = 1
	draft_index_max = length(vessel.weight_dependencies)
	displacement = [w.displacement for w in vessel.weight_dependencies]
	lcb = [w.lcb for w in vessel.weight_dependencies]
	gm = [w.gm for w in vessel.weight_dependencies]
	kg = [w.kg for w in vessel.weight_dependencies]
	metacenter = [w.kmt for w in vessel.weight_dependencies]

	# Get stress limits
	stress_limits = vessel.stress_limits
	shear_min = [s.shear_min for s in stress_limits]
	shear_max = [s.shear_max for s in stress_limits]
	bending_max = [s.bending_max for s in stress_limits]

	# Get frame information
	frame_positions = [f.position for f in vessel.frame_positions]
	# Rasmus: Finds the length of each frame section
	frame_length = diff(frame_positions)  # Assuming uniform frame spacing
	# Rasmus: Find the weight of the ship at each frame section and add a 0 at the beginning of the array
	ship_position_weight = prepend!(diff([f.lightship_weight_cumulative for f in vessel.frame_positions]), 0.0)

	@variable(model, ballast_volume[1:n_ballast_tanks] >= 0)  # Ballast volumes
	@variable(model, cumulative_weight[1:n_positions] >= 0)
	@variable(model, z_min[draft_index_min:draft_index_max], Bin)
	@variable(model, z_max[draft_index_min:draft_index_max], Bin)

	# Force variables
	@variable(model, buoyancy[1:n_positions])
	@variable(model, shear[1:n_positions])
	@variable(model, bending[1:n_positions])

	# Ballast water constraints. Constraint (9)
	@constraint(
		model,
		[t = 1:n_ballast_tanks],
		ballast_volume[t] <= vessel.ballast_tanks[t].max_vol * ρ
	)

	# Draft index constraints. Constraint (3), (4), (5)
	@constraint(model, sum(z_min) == 1)
	@constraint(model, sum(z_max) == 1)
	@constraint(model, [b = draft_index_min:draft_index_max-1], z_min[b] - z_max[b+1] == 0)

	# Calculate ballast tank weights per position. Constraint (10)
	@expression(model, pos_weight_tank[p = 1:n_positions],
    sum(ballast_volume[t] * vessel.ballast_tank_frame_overlap[t, p]
        for t in 1:n_ballast_tanks)
	)

	@expression(model, pos_weight_deadweight[p = 1:n_positions],
		sum(vessel.deadweight_tanks[t].weight * vessel.deadweight_tank_frame_overlap[t, p]
			for t in 1:length(vessel.deadweight_tanks))
	)
	# Modify the total weight calculation to include deadweight tanks. Constraint (1)
	@constraint(model, [p = 1:n_positions],
		cumulative_weight[p] == (p == 1 ? 0 : cumulative_weight[p-1]) +
								pos_weight_cargo[p] +
								pos_weight_tank[p] +
								pos_weight_deadweight[p] +  # Add this line
								ship_position_weight[p]
	)

	# Draft constraints. Constraint (6)
	@constraint(model, dis_min,
		sum(displacement[b] * z_min[b] for b in draft_index_min:draft_index_max) - 
		cumulative_weight[end] <= 0
	)
	# Constraint (7)
	@constraint(model, dis_max,
		sum(displacement[b] * z_max[b] for b in draft_index_min:draft_index_max) - 
		cumulative_weight[end] >= 0
	)

	# Buoyancy calculations using vessel's buoyancy matrix
	# Rasmus: I think this should be constraint (8) but it doesn't look right
	
	@expression(model, buoyancy_interpolated,
		sum(vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_min[b] + 
			vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_max[b] 
			for b in draft_index_min:draft_index_max)./2
	)
	
	#=
	@expression(model, buoyancy_interpolated,
		sum(vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_min[b] 
			for b in draft_index_min:draft_index_max)
	)
	=#

	# # Force relationships
	@constraint(model, buoyancy .== buoyancy_interpolated)

	# Constraint (17) - unsure what is correct
	@constraint(model, shear .== cumulative_weight .- buoyancy)
	#@constraint(model, shear .== cumulative_weight .+ buoyancy)

	# @constraint(model, -100 <= sum(shear[i] * frame_length[i] for i in 1:n_positions-1) <= 100)

	@constraint(model, bending[1] == 0)
	@constraint(model, bending[end] == 0)
	# Constraint (19). 
	# Rasmus: Index for frame i i-1, because of the diff() function to calculate frame_length
	@constraint(model, [i = 2:n_positions],
		bending[i] == bending[i-1] + (shear[i] + shear[i-1]) / 2 * frame_length[i-1]
	)

	stress_positions = [s.position for s in stress_limits]
	# Find the corresponding frame indices for these positions
	stress_frame_indices = [findmin(abs.(frame_positions .- pos))[2] for pos in stress_positions]

	# Stress limits. Constraint (18), (20)
	@constraint(model, [i in 1:length(stress_limits)],
		shear[stress_frame_indices[i]] >= shear_min[i])
	@constraint(model, [i in 1:length(stress_limits)],
		shear[stress_frame_indices[i]] <= shear_max[i])
	@constraint(model, [i in 1:length(stress_limits)],
		bending[stress_frame_indices[i]] <= bending_max[i])

	# # Center of gravity constraints.
	@expression(model, lcg_ballast,
		sum(vessel.ballast_tanks[t].lcg * ballast_volume[t] for t in 1:n_ballast_tanks)
	)
	@expression(model, lcg_deadweight,
    sum(vessel.deadweight_tanks[t].lcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	# Left-hand side of constraint (11) & (12)
	@expression(model, lcg_total,
    	lcg_cargo + lcg_ballast + lcg_deadweight + 
    vessel.lightship_lcg * vessel.lighship_weight
	)
	# Constraint (11)
	@constraint(model, lcg_max,
		lcg_total <= sum(lcb[b] * displacement[b] * z_max[b]
						 for b in draft_index_min:draft_index_max)
	)
	# Constraint (12)
	@constraint(model, lcg_min,
		lcg_total >= sum(lcb[b] * displacement[b] * z_min[b]
						 for b in draft_index_min:draft_index_max)
	)
	# Left-hand side of constraint (13) & (14)
	@expression(model, tcg_ballast,
		sum(vessel.ballast_tanks[t].tcg * ballast_volume[t] for t in 1:n_ballast_tanks)
	)
	@expression(model, tcg_deadweight,
    sum(vessel.deadweight_tanks[t].tcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	@expression(model, tcg_total,
		tcg_cargo + tcg_ballast + tcg_deadweight + 
		vessel.lightship_tcg * vessel.lighship_weight
	)
	# Constraint (13)
	@constraint(model, tcg_max,
		tcg_total <= TCGmax * cumulative_weight[end]
	)
	# Constraint (14)
	@constraint(model, tcg_min,
		tcg_total >= TCGmin * cumulative_weight[end]
	)
	# Some of the left-hand side of constraint (15) & (16)
	# Rasmus: vcg_slope is m_{hat}^{V,T} in paper?
	#@expression(model, vcg_ballast,
	#	sum(vcg_slope[t] * ballast_volume[t] for t in 1:n_ballast_tanks)
	#)
	@variable(model, vcg_ballast >= 0)
	@constraint(model, vcg_ballast == 
		sum(vcg_slope[t] * ballast_volume[t] for t in 1:n_ballast_tanks)
	)
	@expression(model, vcg_deadweight,
    	sum(vessel.deadweight_tanks[t].vcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	# Changed 07/04 
	@expression(model, vcg_total,
		vcg_cargo + vcg_ballast + vcg_deadweight + 
		vessel.lightship_vcg * vessel.lighship_weight #vessel.lightship_vcg
	)

	# # Metacentric height constraint. Constraint (15)
	@constraint(model, metacentric_height,
		sum(gm[b] * displacement[b] * z_min[b] for b in draft_index_min:draft_index_max) <=
		sum(metacenter[b] * displacement[b] * z_min[b] for b in draft_index_min:draft_index_max) - vcg_total
	)

	# KG constraint
	@expression(model, kg_min_total,
		sum(kg[b] * displacement[b] * z_max[b] for b in draft_index_min:draft_index_max)
	)
	# Constraint (16)
	@constraint(model, kg_constraint,
		kg_min_total >= vcg_total
	)

	return model
end


# My slope
function calculate_vcg_slopes_1(vessel::Vessel)
	n_tanks = length(vessel.ballast_tanks)
	slopes = zeros(n_tanks)

	for tank_idx in 1:n_tanks
		# Get tank properties
		tank = vessel.ballast_tanks[tank_idx]
		max_vcg = tank.max_vcg
		max_vcg = tank.max_vcg
		min_vcg = tank.min_vcg
		max_vol = tank.max_vol
		slopes[tank_idx] = round((max_vcg - min_vcg) / max_vol,digits = 2)
		#slopes[tank_idx] = max_vcg
	end
	return slopes
end

function create_model_test(problem::StowageProblem)

	@unpack vessel, slots, cargo = problem

	n_slots = length(slots)
	cargo_types = cargo.cargo_types
	n_positions = length(vessel.frame_positions)
	n_cargo = length(cargo)
	n_deck = length(vessel.decks)
	n_ballast_tanks = length(vessel.ballast_tanks)

	# Rasmus: Fraction of slot within this frame section
	slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)

	# draft_index_min, 
	# draft_index_max,
	# n_positions,

	# Unpack all needed values from common_params
	# @unpack n_ballast_tanks,
	#         ballast_tank_weight_max,
	#         ballast_tank_proportions,
	#         ship_position_weight,
	#         ship_buoyancy_weight,
	#         ship_buoyancy_position,
	#         force_index,
	#         shear_min,
	#         shear_max,
	#         bending_max = problem.vessel

	# Unpack needed values from solver struct m
	# @unpack parameter, ship = m
	# @unpack CSC, D, draft_index_full, TCGmax, TCGmin, vcgnum, lcgnum, B = parameter
	# @unpack regulartank, hydro, deck = ship
	# @unpack lcg, tcg = regulartank
	# @unpack displacement, lcb, gm_min, metacenter, kg_min = hydro

	ρ = 1.025

	#model = Model(() -> HiGHS.Optimizer())
	#set_time_limit_sec(model, 60 * 0.5)
	# model = Model(() -> Gurobi.Optimizer(GRB_ENV[]))
	model = Model(Gurobi.Optimizer)
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)
	# model = Model(HiGHS.Optimizer)

	# Weight variables
	@variable(model, weight[1:n_slots] >= 0)   # Weight at each slot

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

	for c in eachindex(cargo)
		for s in eachindex(slots)
			# If assignment is incompatible, force x[i,j] to be 0
			if !check_refrigeration_compatibility(cargo[c], slots[s])
				fix(cs[c, s], 0; force = true)
			end
		end
	end

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

	# Defining cargo_slack variable. Constraint (23)
	@constraint(
		model,
		[c = 1:n_cargo],
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

	# Rasmus: expression defines something (pos_weight_cargo) without it being a variable or constraint
	# Rasmus: pos_weight_cargo is used as input later
	# Rasmus: This calculates the weight of cargo in each frame. Constraint (28)
	@expression(model, pos_weight_cargo[p = 1:n_positions],
		# Cargo weights using precalculated proportions
		sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots)
	)

	# # # Deck weight limit. Constraint (29)
	@constraint(
		model,
		[d = 1:n_deck],
		sum(weight[s] for s in [x.id for x in filter(x -> x.deck_id == d, slots)]) <= vessel.decks[d].weight_limit
	)

	# Rasmus: lcg, vcg, and tcg for cargo. Constraint (30), (31), (32)
	@expression(model, lcg_cargo, sum(weight[s] * slots[s].lcg for s ∈ 1:n_slots))
	@expression(model, tcg_cargo, sum(weight[s] * slots[s].tcg for s ∈ 1:n_slots))
	# Rasmus: Don't understand why is not using weight[s]
	@expression(model, vcg_cargo, sum(cs[c, s] * cargo[c].weight * (slots[s].vcg + get_height(cargo[c]) / 2) for s ∈ 1:n_slots, c ∈ 1:n_cargo))

	# Rasmus: y_t is defined in this part of the model
	add_stability_test!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)

	ballast_volume = model[:ballast_volume]

	@objective(model, Min,
		sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
		-
		sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
	)

	return model
end

function add_stability_test!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)
	vcg_slope = calculate_vcg_slopes_1(vessel)

	n_positions = length(vessel.frame_positions)
	n_ballast_tanks = length(vessel.ballast_tanks)

	# Constants
	ρ = 1.025  # Water density
	TCGmax = 0.001  # TODO: Get from vessel
	TCGmin = -0.001  # TODO: Get from vessel

	# Extract weight dependencies
	draft_index_min = 1
	draft_index_max = length(vessel.weight_dependencies)
	displacement = [w.displacement for w in vessel.weight_dependencies]
	lcb = [w.lcb for w in vessel.weight_dependencies]
	gm = [w.gm for w in vessel.weight_dependencies]
	kg = [w.kg for w in vessel.weight_dependencies]
	metacenter = [w.kmt for w in vessel.weight_dependencies]

	# Get stress limits
	stress_limits = vessel.stress_limits
	shear_min = [s.shear_min for s in stress_limits]
	shear_max = [s.shear_max for s in stress_limits]
	bending_max = [s.bending_max for s in stress_limits]

	# Get frame information
	frame_positions = [f.position for f in vessel.frame_positions]
	# Rasmus: Finds the length of each frame section
	frame_length = diff(frame_positions)  # Assuming uniform frame spacing
	# Rasmus: Find the weight of the ship at each frame section and add a 0 at the beginning of the array
	ship_position_weight = prepend!(diff([f.lightship_weight_cumulative for f in vessel.frame_positions]), 0.0)

	@variable(model, ballast_volume[1:n_ballast_tanks] >= 0)  # Ballast volumes
	@variable(model, cumulative_weight[1:n_positions] >= 0)
	@variable(model, z_min[draft_index_min:draft_index_max], Bin)
	@variable(model, z_max[draft_index_min:draft_index_max], Bin)

	# Force variables
	@variable(model, buoyancy[1:n_positions])
	@variable(model, shear[1:n_positions])
	@variable(model, bending[1:n_positions])

	# Ballast water constraints. Constraint (9)
	@constraint(
		model,
		[t = 1:n_ballast_tanks],
		ballast_volume[t] <= vessel.ballast_tanks[t].max_vol * ρ
	)

	# Draft index constraints. Constraint (3), (4), (5)
	@constraint(model, sum(z_min) == 1)
	@constraint(model, sum(z_max) == 1)
	@constraint(model, [b = draft_index_min:draft_index_max-1], z_min[b] - z_max[b+1] == 0)

	# Calculate ballast tank weights per position. Constraint (10)
	@expression(model, pos_weight_tank[p = 1:n_positions],
    sum(ballast_volume[t] * vessel.ballast_tank_frame_overlap[t, p]
        for t in 1:n_ballast_tanks)
	)

	@expression(model, pos_weight_deadweight[p = 1:n_positions],
		sum(vessel.deadweight_tanks[t].weight * vessel.deadweight_tank_frame_overlap[t, p]
			for t in 1:length(vessel.deadweight_tanks))
	)
	# Modify the total weight calculation to include deadweight tanks. Constraint (1)
	@constraint(model, [p = 1:n_positions],
		cumulative_weight[p] == (p == 1 ? 0 : cumulative_weight[p-1]) +
								pos_weight_cargo[p] +
								pos_weight_tank[p] +
								pos_weight_deadweight[p] +  # Add this line
								ship_position_weight[p]
	)

	# Draft constraints. Constraint (6)
	@constraint(model, dis_min,
		sum(displacement[b] * z_min[b] for b in draft_index_min:draft_index_max) - 
		cumulative_weight[end] <= 0
	)
	# Constraint (7)
	@constraint(model, dis_max,
		sum(displacement[b] * z_max[b] for b in draft_index_min:draft_index_max) - 
		cumulative_weight[end] >= 0
	)

	# Buoyancy calculations using vessel's buoyancy matrix
	# Rasmus: I think this should be constraint (8) but it doesn't look right
	
	@expression(model, buoyancy_interpolated,
		sum(vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_min[b] + 
			vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_max[b] 
			for b in draft_index_min:draft_index_max)./2
	)
	
	#=
	@expression(model, buoyancy_interpolated,
		sum(vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_min[b] 
			for b in draft_index_min:draft_index_max)
	)
	=#

	# # Force relationships
	@constraint(model, buoyancy .== buoyancy_interpolated)

	# Constraint (17) - unsure what is correct
	@constraint(model, shear .== cumulative_weight .- buoyancy)
	#@constraint(model, shear .== cumulative_weight .+ buoyancy)

	# @constraint(model, -100 <= sum(shear[i] * frame_length[i] for i in 1:n_positions-1) <= 100)

	@constraint(model, bending[1] == 0)
	@constraint(model, bending[end] == 0)
	# Constraint (19). 
	# Rasmus: Index for frame i i-1, because of the diff() function to calculate frame_length
	@constraint(model, [i = 2:n_positions],
		bending[i] == bending[i-1] + (shear[i] + shear[i-1]) / 2 * frame_length[i-1]
	)

	stress_positions = [s.position for s in stress_limits]
	# Find the corresponding frame indices for these positions
	stress_frame_indices = [findmin(abs.(frame_positions .- pos))[2] for pos in stress_positions]

	# Stress limits. Constraint (18), (20)
	@constraint(model, [i in 1:length(stress_limits)],
		shear[stress_frame_indices[i]] >= shear_min[i])
	@constraint(model, [i in 1:length(stress_limits)],
		shear[stress_frame_indices[i]] <= shear_max[i])
	@constraint(model, [i in 1:length(stress_limits)],
		bending[stress_frame_indices[i]] <= bending_max[i])

	# # Center of gravity constraints.
	@expression(model, lcg_ballast,
		sum(vessel.ballast_tanks[t].lcg * ballast_volume[t] for t in 1:n_ballast_tanks)
	)
	@expression(model, lcg_deadweight,
    sum(vessel.deadweight_tanks[t].lcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	# Left-hand side of constraint (11) & (12)
	@expression(model, lcg_total,
    	lcg_cargo + lcg_ballast + lcg_deadweight + 
    vessel.lightship_lcg * vessel.lighship_weight
	)
	# Constraint (11)
	@constraint(model, lcg_max,
		lcg_total <= sum(lcb[b] * displacement[b] * z_max[b]
						 for b in draft_index_min:draft_index_max)
	)
	# Constraint (12)
	@constraint(model, lcg_min,
		lcg_total >= sum(lcb[b] * displacement[b] * z_min[b]
						 for b in draft_index_min:draft_index_max)
	)
	# Left-hand side of constraint (13) & (14)
	@expression(model, tcg_ballast,
		sum(vessel.ballast_tanks[t].tcg * ballast_volume[t] for t in 1:n_ballast_tanks)
	)
	@expression(model, tcg_deadweight,
    sum(vessel.deadweight_tanks[t].tcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	@expression(model, tcg_total,
		tcg_cargo + tcg_ballast + tcg_deadweight + 
		vessel.lightship_tcg * vessel.lighship_weight
	)
	# Constraint (13)
	@constraint(model, tcg_max,
		tcg_total <= TCGmax * cumulative_weight[end]
	)
	# Constraint (14)
	@constraint(model, tcg_min,
		tcg_total >= TCGmin * cumulative_weight[end]
	)
	# Some of the left-hand side of constraint (15) & (16)
	# Rasmus: vcg_slope is m_{hat}^{V,T} in paper?
	#@expression(model, vcg_ballast,
	#	sum(vcg_slope[t] * ballast_volume[t] for t in 1:n_ballast_tanks)
	#)
	@variable(model, vcg_ballast >= 0)
	@constraint(model, vcg_ballast == 
		sum(vcg_slope[t] * ballast_volume[t] for t in 1:n_ballast_tanks)
	)
	@expression(model, vcg_deadweight,
    	sum(vessel.deadweight_tanks[t].vcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	# Changed 07/04 
	@expression(model, vcg_total,
		vcg_cargo + vcg_ballast + vcg_deadweight + 
		vessel.lightship_vcg * vessel.lighship_weight #vessel.lightship_vcg
	)

	# # Metacentric height constraint. Constraint (15)
	@constraint(model, metacentric_height,
		sum(gm[b] * displacement[b] * z_min[b] for b in draft_index_min:draft_index_max) <=
		sum(metacenter[b] * displacement[b] * z_min[b] for b in draft_index_min:draft_index_max) - vcg_total
	)

	# KG constraint
	@expression(model, kg_min_total,
		sum(kg[b] * displacement[b] * z_max[b] for b in draft_index_min:draft_index_max)
	)
	# Constraint (16)
	@constraint(model, kg_constraint,
		kg_min_total >= vcg_total
	)

	return model
end