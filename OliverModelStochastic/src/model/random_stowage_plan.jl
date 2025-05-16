# Creates random stowage plan and minimizes shifts and ballast-water until feasible

# Generates random stowage plan
function random_stowage_plan(CargoC::CargoCollection, slots)
    # Function creates a list of slot ids which are equal to the cargo type id input "t"
    valid_slots(t::Int) = [slot.id for slot in filter(x -> x.cargo_type_id == t, slots)]
    n_slots = length(slots)
    n_cargo = length(CargoC.items)
    # stowage plan
    cs = zeros(Bool, n_cargo, n_slots)
    not_stowaged = []
    # Valid slots with respect to slots already taken
    slot_overlap = zeros(Int, n_slots, n_slots)
    n_placed = 0
    for i in 1:n_cargo
        placed = false
        val_slots_type = valid_slots(CargoC[i].cargo_type_id)
        set_slots = []
        for j in val_slots_type
            if sum(slot_overlap[j,:]) < 0.5
                push!(set_slots, j)
            end
        end
        if length(set_slots) == 0 # not any feasible spots left
            push!(not_stowaged, CargoC[i])
            continue
        else
            slot_id = rand(set_slots)
            cs[CargoC[i].id, slot_id] = 1
            slot_overlap[:,slot_id] = copy(Int.(slots.overlap_matrix[slot_id,:]))
            slot_overlap[slot_id,:] = copy(Int.(slots.overlap_matrix[slot_id,:]))
            slot_overlap[slot_id,slot_id] = 1
            placed = true
            n_placed += 1
        end
        if !placed
            push!(not_stowaged, CargoC[i])
        end
    end
    #println("Placed: ", n_placed)
    return cs, not_stowaged
end

# Generate random cargocollection given how much cargo is wanted.
function random_cargocollection(ntypes::Vector{Int}, shuffle::Bool = true)
    cargoes = Cargo[]
    id = 0
    for i in 1:4
        for j in 1:ntypes[i]
            #cargo_types = [Car, Truck, Machine, Secu]
            if i == 1 # truck
                cargo_type = Truck
            elseif i == 2 # cargo
                cargo_type = Car
            elseif i == 3 # machine
                cargo_type = Machine
            else # secu
                cargo_type = Secu
            end
            type_info = CARGO_TYPES[string(cargo_type)]
    
            # Weight based on cargo type and volume
            volume = type_info.length * type_info.width * type_info.height
    
            weight = if cargo_type == Car
            rand(1.5:0.1:3.0)
            elseif cargo_type == Truck
                rand(5.0:0.5:40.0)
            else
                volume * rand(0.2:0.01:0.5) # density factor
            end
            loading_port = 1
            discharge_port = 2
            id += 1
            push!(cargoes,Cargo(
                id=id,
                cargo_type_id=type_info.id,  # Use the id from type_info
                weight=weight,
                loading_port=loading_port,
                discharge_port=discharge_port,
                priority=1,#rand(1:3),
                requires_lashing=false,#rand(Bool),
                requires_ventilation=false,#rand(Bool),
                hazardous=rand(0:3), # doesn't matter
                refers=false#rand(Bool)
            ))
        end
    end
    if shuffle 
        return CargoCollection(shuffle!(cargoes))
    else
        return CargoCollection(cargoes)
    end
end

# sort cargo collection by type
function sort_cargocollection(cargoc::CargoCollection, order::Vector{Int64})
    cargoes = [cargoc.items[i] for i in 1:length(cargoc.items)]
    # Sort by cargo type - change order if wanted
    sorted_cargo = filter(x -> x.cargo_type_id == order[1], cargoes)
    append!(sorted_cargo,filter(x -> x.cargo_type_id == order[2], cargoes))
    append!(sorted_cargo,filter(x -> x.cargo_type_id == order[3], cargoes))
    append!(sorted_cargo,filter(x -> x.cargo_type_id == order[4], cargoes))
    return CargoCollection(sorted_cargo)
end

# Minimizes ballast water and shifts. 
# Only shifts if it cannot stabilize the ship otherwise, 
# i.e. if random stowage plan is not feasible
function create_random_stowageplan_model(cs_old,not_stowed,cargo,vessel,slots, new_cargo_allowed::Bool = false)

    n_slots = length(slots)
	cargo_types = cargo.cargo_types
	n_positions = length(vessel.frame_positions)
	n_cargo = length(cargo)
	n_deck = length(vessel.decks)
	n_ballast_tanks = length(vessel.ballast_tanks)
	# Fraction of slot within this frame section
	slots_to_frame = calculate_slot_frame_overlap_matrix(slots, vessel.frame_positions)
    # Water density
	ρ = 1.025
    # Model
	model = Model(Gurobi.Optimizer)
	# number should match number of cores used at HPC
	set_optimizer_attribute(model, "Threads", 4)

    # Weight variables
	@variable(model, weight[1:n_slots] >= 0)   # Weight at each slot
	# Cargo assignment
	# Rasmus: cs is a binary matrix, 1 if cargo c is assigned to slot s
	@variable(model, cs[1:n_cargo, 1:n_slots], Bin)  # Assignment variables
	@variable(model, cargo_slack[1:n_cargo], Bin) # 1 if cargo is assigned a slot
    @variable(model, y[1:n_cargo], Bin) # 1 if cargo i is moved from orignal stowage plan

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
        sum(cs[c, s] for s ∈ 1:n_slots) == cargo_slack[c])
    #   Rasmus: Function creates a list of slot ids which are not equal to the cargo type id input "t"
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
		sum(weight[s] * slots_to_frame[s, p] for s ∈ 1:n_slots))

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

    # If we do not allow new cargo to be stowed
    if !new_cargo_allowed
        @constraint(model, [i in [cargo.id for cargo in not_stowed]], cargo_slack[i] == 0)
    end

    # Penalty for moving cargo
	M = 100000 # Should be determined more precisely at some point
    #M = sum(cost)+1

    @constraint(model, [c = 1:n_cargo, s = 1:n_slots], 
        y[c] >= cs[c,s] - cs_old[c,s]) # cargo is moved
    @constraint(model, [c = 1:n_cargo, s = 1:n_slots], 
        y[c] >= cs_old[c,s] - cs[c,s]) # cargo is moved

    @objective(model, Min,
		sum(ballast_volume[t] for t ∈ 1:n_ballast_tanks)
		-
		sum(cost[c] * cargo_slack[c] for c ∈ 1:n_cargo)
        + M*sum(y)  # penalty for moving cargo
	)
	return model
end


