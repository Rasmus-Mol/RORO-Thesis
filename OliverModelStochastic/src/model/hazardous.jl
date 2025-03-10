
function add_hazardous_cargo_constraints!(model, vessel, cargos, m)

    cs = model[:cs]
    n_cargo = length(cargos)
    n_slots = length(vessel.slots)
    
    @expression(model, hazardous_cargo, 
        sum(cargos[c].haz_class * cs[c, s] 
        for c = 1:n_cargo, 
            s = 1:n_slots))

    MOI.set(
        model,
        MOI.LazyConstraintCallback(),
        (cb_data) -> hazardous_lazy_constraints(cb_data, vessel, model, cs, cargos, m)
    )

    MOI.set(
        model,
        MOI.HeuristicCallback(),
        (cb_data) -> hazardous_heuristic(cb_data, vessel, model, cs, cargos, m)
    )

end

function calc_slot_dist(s1::Slot, s2::Slot)
    lcg_dist = abs(s1.lcg - s2.lcg)
    tcg_dist = abs(s1.tcg - s2.tcg)
    return (lcg_dist, tcg_dist)
end

function get_corners(s::Slot)
    # Assuming width and length properties exist
    half_width = s.width / 2
    half_length = s.length / 2
    
    return [
        (s.lcg - half_length, s.tcg - half_width),  # Bottom left
        (s.lcg - half_length, s.tcg + half_width),  # Top left
        (s.lcg + half_length, s.tcg - half_width),  # Bottom right
        (s.lcg + half_length, s.tcg + half_width)   # Top right
    ]
end

function calc_euclidean_dist(s1::Slot, s2::Slot)
    corners1 = get_corners(s1)
    corners2 = get_corners(s2)
    
    min_dist = Inf
    for c1 in corners1
        for c2 in corners2
            dist = sqrt((c1[1] - c2[1])^2 + (c1[2] - c2[2])^2)
            min_dist = min(min_dist, dist)
        end
    end
    return min_dist
end


function calc_all_slot_distances(slots::Vector{Slot})
    n = length(slots)
    lcg_distances = zeros(Float64, n, n)
    tcg_distances = zeros(Float64, n, n)
    
    for i in 1:n
        for j in (i+1):n
            lcg_dist, tcg_dist = calc_slot_dist(slots[i], slots[j])
            lcg_distances[i,j] = lcg_dist
            lcg_distances[j,i] = lcg_dist
            tcg_distances[i,j] = tcg_dist
            tcg_distances[j,i] = tcg_dist
        end
    end
    
    return (lcg_distances, tcg_distances)
end

function check_violation(m, danger_rule, lcgdist, tcgdist, deck_diff, same_deck, on_deck)
    if danger_rule == 2
        return same_deck && all([tcgdist, lcgdist] .< m.parameter.mindist[1:2])
    elseif danger_rule == 3
        if !on_deck
            return same_deck || all([tcgdist, lcgdist] .< m.parameter.mindist[4])
        else
            return same_deck && all([tcgdist, lcgdist] .< m.parameter.mindist[3])
        end
    elseif danger_rule >= 4  # danger_rule is 4 or 5
        if !on_deck
            return deck_diff < 2 || lcgdist < m.parameter.mindist[5]
        else
            return same_deck && lcgdist < m.parameter.mindist[5]
        end
    end
    return false
end

function hazardous_lazy_constraints(cb_data, vessel, model, cs, cargos, m)
    # Get placement indexes
    cs_val = callback_value.(Ref(cb_data), cs)
    placed_cargo = findall(cs_val .> 0.5)
    
    hazardous_cargo_ids = [cargo.id for cargo in cargos if cargo.haz_class != 18]

    # Filter for dangerous cargo
    placed_cargo_hazardous = [Tuple(idx) for idx in placed_cargo if idx[1] in hazardous_cargo_ids]

    n_lazy_constraints = 0

    for i in eachindex(placed_cargo_hazardous[1:end])
        cargo1_id, slot1_id = placed_cargo_hazardous[i]
        cargo1_class = cargos[cargo1_id].haz_class
        slot1 = vessel.slots[slot1_id].slot
        
        for j in eachindex(placed_cargo_hazardous[1:end])
            if i == j
                continue
            end
            cargo2_id, slot2_id = placed_cargo_hazardous[j]
            cargo2_class = cargos[cargo2_id].haz_class
            slot2 = vessel.slots[slot2_id].slot
            
            # Check danger rule
            danger_rule = m.parameter.segtable[cargo1_class, cargo2_class]
            if danger_rule < 2
                continue
            end
            
            # Calculate deck properties
            deck1 = vessel.slots[slot1_id].deck_id
            deck2 = vessel.slots[slot2_id].deck_id
            deck_diff = abs(deck1 - deck2)
            same_deck = deck_diff == 0
            on_deck = (slot1.on_deck) && (slot2.on_deck)

            dist = calc_euclidean_dist(slot1, slot2)
            dist_lcg = dist
            dist_tcg = dist

            # dist_lcg, dist_tcg = calc_slot_dist(slot1, slot2) # TODO Precalculate

            # dist_lcg = dists_lcg[slot1_id, slot2_id]
            # dist_tcg = dists_tcg[slot1_id, slot2_id]
            
            # Skip if one cargo on deck
            # (on_deck == 1) && continue
            
            # Skip if satisfies strongest rule
            if deck_diff >= 2 && dist_lcg >= m.parameter.mindist[5]
                continue
            end            
     
            # Check for violations
            if check_violation(m, danger_rule, dist_lcg, dist_tcg, deck_diff, same_deck, on_deck)
                MOI.submit(
                    model,
                    MOI.LazyConstraint(cb_data),
                    @build_constraint(cs[cargo1_id, slot1_id] + cs[cargo2_id, slot2_id] <= 1)
                )
                n_lazy_constraints += 1
            end
        end
    end
end

mutable struct HeuristicStatus
    count::Int
    status::Any
end

const HEURISTIC_STATUS = HeuristicStatus(0, "")


function hazardous_heuristic(cb_data, vessel, model, cs, cargos, m)
    if HEURISTIC_STATUS.status == MOI.HEURISTIC_SOLUTION_ACCEPTED || HEURISTIC_STATUS.count >= 10
        return
    else
        println("Heuristic called $(HEURISTIC_STATUS.count) times")
        status = callback_node_status(cb_data, model)
    # if status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        # Get current solution values
        println(status)
        cs_val = callback_value.(Ref(cb_data), cs)
        # Create copy of solution
        new_sol = copy(cs_val)
        
        # Find hazardous cargo indices
        hazardous_indices = findall(x->x.haz_class != 18, cargos)
        
        # Set solution values to 0 for hazardous cargo
        for i in hazardous_indices
            for j in 1:size(cs)[2]
                new_sol[i,j] = 0
            end
        end
        
        # Submit modified solution - flatten variables and values
        status = MOI.submit(
            model,
            MOI.HeuristicSolution(cb_data),
            vec(cs),
            vec(new_sol)
        )
        println("Submitted heuristic solution with status: ", status)
        HEURISTIC_STATUS.count += 1  # Increment counter
        HEURISTIC_STATUS.status = status
        return
    end
end

function get_heuristic_call_count()
    return HEURISTIC_COUNTER.count
end