# Rasmus: Creates an array with slope vcg for all tanks
function calculate_vcg_slopes(vessel::Vessel)
	n_tanks = length(vessel.ballast_tanks)
	slopes = zeros(n_tanks)
	for tank_idx in 1:n_tanks
		# Get tank properties
		tank = vessel.ballast_tanks[tank_idx]
		max_vcg = tank.max_vcg
		min_vcg = tank.min_vcg
		max_vol = tank.max_vol
		slopes[tank_idx] = (max_vcg - min_vcg) / max_vol
	end
	return slopes
end

# Rasmus: Add stability constraints to the model
# pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo are scenario dependent
function add_stability_stochastic!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo, scenarios)
	vcg_slope = calculate_vcg_slopes(vessel)

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

	@variable(model, ballast_volume[1:n_ballast_tanks,1:scenarios] >= 0)  # Ballast volumes
	@variable(model, cumulative_weight[1:n_positions,1:scenarios] >= 0)
	# Unsure if z should be two-stage or first-stage
	@variable(model, z_min[draft_index_min:draft_index_max,1:scenarios], Bin)
	@variable(model, z_max[draft_index_min:draft_index_max,1:scenarios], Bin)

	# Force variables
	# Unsure if they should be two-stage or first-stage
	@variable(model, buoyancy[1:n_positions,1:scenarios])
	@variable(model, shear[1:n_positions,1:scenarios])
	@variable(model, bending[1:n_positions,1:scenarios])

	# Ballast water constraints. Constraint (9)
	@constraint(
		model,
		[t = 1:n_ballast_tanks,sc=1:scenarios],
		ballast_volume[t,sc] <= vessel.ballast_tanks[t].max_vol * ρ
	)

	# Draft index constraints. Constraint (3), (4), (5)
	@constraint(model, [sc = 1:scenarios], sum(z_min[i,sc] for i in draft_index_min:draft_index_max) == 1)
	@constraint(model, [sc = 1:scenarios], sum(z_max[i,sc] for i in draft_index_min:draft_index_max) == 1)
	@constraint(model, [b = draft_index_min:draft_index_max-1,sc=1:scenarios], z_min[b,sc] - z_max[b+1,sc] == 0)

	# Calculate ballast tank weights per position. Constraint (10)
	@expression(model, pos_weight_tank[p = 1:n_positions, sc = 1:scenarios],
    sum(ballast_volume[t,sc] * vessel.ballast_tank_frame_overlap[t, p]
        for t in 1:n_ballast_tanks)
	)
	@expression(model, pos_weight_deadweight[p = 1:n_positions],
		sum(vessel.deadweight_tanks[t].weight * vessel.deadweight_tank_frame_overlap[t, p]
			for t in 1:length(vessel.deadweight_tanks))
	)
	# Modify the total weight calculation to include deadweight tanks. Constraint (1)
	@constraint(model, [p = 1:n_positions, sc = 1:scenarios],
		cumulative_weight[p,sc] == (p == 1 ? 0 : cumulative_weight[p-1,sc]) +
								pos_weight_cargo[p,sc] +
								pos_weight_tank[p,sc] +
								pos_weight_deadweight[p] +  # Add this line
								ship_position_weight[p]
	)
	# Draft constraints. Constraint (6)
	@constraint(model, dis_min[sc = 1:scenarios],
		sum(displacement[b] * z_min[b,sc] for b in draft_index_min:draft_index_max) - 
		cumulative_weight[end,sc] <= 0
	)
	# Constraint (7)
	@constraint(model, dis_max[sc = 1:scenarios],
		sum(displacement[b] * z_max[b,sc] for b in draft_index_min:draft_index_max) - 
		cumulative_weight[end,sc] >= 0
	)
	# Buoyancy calculations using vessel's buoyancy matrix
	# Rasmus: I think this should be constraint (8) but it doesn't look right
	@expression(model, buoyancy_interpolated[sc = 1:scenarios],
		sum(vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_min[b,sc] +
			vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_max[b,sc] 
			for b in draft_index_min:draft_index_max)
	)
	# # Force relationships
	@constraint(model, [sc = 1:scenarios], buoyancy[:,sc] .== buoyancy_interpolated[sc])

	# Constraint (17) - Unsure what is correct
	@constraint(model, shear .== cumulative_weight .- buoyancy)
	#@constraint(model, shear .== cumulative_weight .+ buoyancy)

	# @constraint(model, -100 <= sum(shear[i] * frame_length[i] for i in 1:n_positions-1) <= 100)

	@constraint(model,[sc = 1:scenarios], bending[1,sc] == 0)
	@constraint(model,[sc = 1:scenarios], bending[end,sc] == 0)
	# Constraint (19). 
	# Rasmus: Index for frame i i-1, because of the diff() function to calculate frame_length
	@constraint(model, [i = 2:n_positions, sc = 1:scenarios],
		bending[i,sc] == bending[i-1,sc] + (shear[i,sc] + shear[i-1,sc]) / 2 * frame_length[i-1]
	)
	stress_positions = [s.position for s in stress_limits]
	# Find the corresponding frame indices for these positions
	stress_frame_indices = [findmin(abs.(frame_positions .- pos))[2] for pos in stress_positions]

	# Stress limits. Constraint (18), (20)
	@constraint(model, [i in 1:length(stress_limits),sc = 1:scenarios],
		shear[stress_frame_indices[i],sc] >= shear_min[i])
	@constraint(model, [i in 1:length(stress_limits),sc = 1:scenarios],
		shear[stress_frame_indices[i],sc] <= shear_max[i])
	@constraint(model, [i in 1:length(stress_limits),sc = 1:scenarios],
		bending[stress_frame_indices[i],sc] <= bending_max[i])

	# # Center of gravity constraints.
	@expression(model, lcg_ballast[sc = 1:scenarios],
		sum(vessel.ballast_tanks[t].lcg * ballast_volume[t,sc] for t in 1:n_ballast_tanks)
	)
	@expression(model, lcg_deadweight,
    sum(vessel.deadweight_tanks[t].lcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	# Left-hand side of constraint (11) & (12)
	@expression(model, lcg_total[sc = 1:scenarios],
    	lcg_cargo[sc] + lcg_ballast[sc] + lcg_deadweight + 
    vessel.lightship_lcg * vessel.lighship_weight
	)
	# Constraint (11)
	@constraint(model, lcg_max[sc = 1:scenarios],
		lcg_total[sc] <= sum(lcb[b] * displacement[b] * z_max[b,sc]
						 for b in draft_index_min:draft_index_max)
	)
	# Constraint (12)
	@constraint(model, lcg_min[sc = 1:scenarios],
		lcg_total[sc] >= sum(lcb[b] * displacement[b] * z_min[b,sc]
						 for b in draft_index_min:draft_index_max)
	)
	# Left-hand side of constraint (13) & (14)
	@expression(model, tcg_ballast[sc = 1:scenarios],
		sum(vessel.ballast_tanks[t].tcg * ballast_volume[t,sc] for t in 1:n_ballast_tanks)
	)
	@expression(model, tcg_deadweight,
    sum(vessel.deadweight_tanks[t].tcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	@expression(model, tcg_total[sc = 1:scenarios],
		tcg_cargo[sc] + tcg_ballast[sc] + tcg_deadweight + 
		vessel.lightship_tcg * vessel.lighship_weight
	)
	# Constraint (13)
	@constraint(model, tcg_max[sc = 1:scenarios],
		tcg_total[sc] <= TCGmax * cumulative_weight[end]
	)
	# Constraint (14)
	@constraint(model, tcg_min[sc = 1:scenarios],
		tcg_total[sc] >= TCGmin * cumulative_weight[end]
	)
	# Some of the left-hand side of constraint (15) & (16)
	# Rasmus: vcg_slope relates to m_{hat}^{V,T} in the paper
	@expression(model, vcg_ballast[sc = 1:scenarios],
		sum(vcg_slope[t] * ballast_volume[t,sc] for t in 1:n_ballast_tanks)
	)
	@expression(model, vcg_deadweight,
    	sum(vessel.deadweight_tanks[t].vcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	# Changed 07/04
	@expression(model, vcg_total[sc = 1:scenarios],
		vcg_cargo[sc] + vcg_ballast[sc] + vcg_deadweight + 
		vessel.lightship_vcg * vessel.lighship_weight #vessel.lightship_vcg
	)

	# # Metacentric height constraint. Constraint (15) - NB DIFFERENT THAN IN ARTICLE
	@constraint(model, metacentric_height[sc = 1:scenarios],
		sum(gm[b] * displacement[b] * z_min[b,sc] for b in draft_index_min:draft_index_max) <=
		sum(metacenter[b] * displacement[b] * z_min[b,sc] for b in draft_index_min:draft_index_max) - vcg_total[sc]
	)

	# KG constraint - NB DIFFERENT THAN IN ARTICLE
	@expression(model, kg_min_total[sc = 1:scenarios],
		sum(kg[b] * displacement[b] * z_max[b,sc] for b in draft_index_min:draft_index_max)
	)
	# Constraint (16)
	@constraint(model, kg_constraint[sc = 1:scenarios],
		kg_min_total[sc] >= vcg_total[sc]
	)
	return model
end


# Rasmus: Add stability constraints to the model
# pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo are scenario dependent
function add_stability_stochastic_slack_all!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo, scenarios, M)
	vcg_slope = calculate_vcg_slopes(vessel)

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

	@variable(model, ballast_volume[1:n_ballast_tanks,1:scenarios] >= 0)  # Ballast volumes
	@variable(model, cumulative_weight[1:n_positions,1:scenarios] >= 0)
	@variable(model, z_min[draft_index_min:draft_index_max,1:scenarios], Bin)
	@variable(model, z_max[draft_index_min:draft_index_max,1:scenarios], Bin)

	# Force variables
	@variable(model, buoyancy[1:n_positions,1:scenarios])
	@variable(model, shear[1:n_positions,1:scenarios])
	@variable(model, bending[1:n_positions,1:scenarios])

	# Slack variables:
	# Center of gravity
	@variable(model, slack_Vmax[1:scenarios] >= 0)
	@variable(model, slack_Vmin[1:scenarios] >= 0)
	@variable(model, slack_Tmin[1:scenarios] >= 0)
	@variable(model, slack_Tmax[1:scenarios] >= 0)
	@variable(model, slack_Lmin[1:scenarios] >= 0)
	@variable(model, slack_Lmax[1:scenarios] >= 0)
	# Stress and Bending
	#@variable(model, slack_shear1[1:n_positions,1:scenarios]>=0)
	#@variable(model, slack_shear2[1:n_positions,1:scenarios]>=0)
	stress_limits = vessel.stress_limits
	@variable(model, slack_shearMin[1:length(stress_limits),1:scenarios] >= 0)
	@variable(model, slack_shearMax[1:length(stress_limits),1:scenarios] >= 0)
	@variable(model, slack_bendingMax[1:length(stress_limits),1:scenarios] >= 0)
	# ballast tanks
	@variable(model, slack_ballast_tanks[1:n_ballast_tanks,1:scenarios] >= 0)

	# Ballast water constraints. Constraint (9)
	@constraint(
		model,
		[t = 1:n_ballast_tanks,sc=1:scenarios],
		ballast_volume[t,sc] <= vessel.ballast_tanks[t].max_vol * ρ + slack_ballast_tanks[t,sc]
	)

	# Draft index constraints. Constraint (3), (4), (5)
	@constraint(model, [sc = 1:scenarios], sum(z_min[i,sc] for i in draft_index_min:draft_index_max) == 1)
	@constraint(model, [sc = 1:scenarios], sum(z_max[i,sc] for i in draft_index_min:draft_index_max) == 1)
	@constraint(model, [b = draft_index_min:draft_index_max-1,sc=1:scenarios], z_min[b,sc] - z_max[b+1,sc] == 0)

	# Calculate ballast tank weights per position. Constraint (10)
	@expression(model, pos_weight_tank[p = 1:n_positions, sc = 1:scenarios],
    sum(ballast_volume[t,sc] * vessel.ballast_tank_frame_overlap[t, p]
        for t in 1:n_ballast_tanks)
	)
	@expression(model, pos_weight_deadweight[p = 1:n_positions],
		sum(vessel.deadweight_tanks[t].weight * vessel.deadweight_tank_frame_overlap[t, p]
			for t in 1:length(vessel.deadweight_tanks))
	)
	# Modify the total weight calculation to include deadweight tanks. Constraint (1)
	@constraint(model, [p = 1:n_positions, sc = 1:scenarios],
		cumulative_weight[p,sc] == (p == 1 ? 0 : cumulative_weight[p-1,sc]) +
								pos_weight_cargo[p,sc] +
								pos_weight_tank[p,sc] +
								pos_weight_deadweight[p] +  # Add this line
								ship_position_weight[p]
	)
	# Draft constraints. Constraint (6)
	@constraint(model, dis_min[sc = 1:scenarios],
		sum(displacement[b] * z_min[b,sc] for b in draft_index_min:draft_index_max) - 
		cumulative_weight[end,sc] <= 0
	)
	# Constraint (7)
	@constraint(model, dis_max[sc = 1:scenarios],
		sum(displacement[b] * z_max[b,sc] for b in draft_index_min:draft_index_max) - 
		cumulative_weight[end,sc] >= 0
	)
	# Buoyancy calculations using vessel's buoyancy matrix
	# Rasmus: I think this should be constraint (8) but it doesn't look right
	@expression(model, buoyancy_interpolated[sc = 1:scenarios],
		sum(vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_min[b,sc] +
			vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_max[b,sc]
			for b in draft_index_min:draft_index_max)
	)
	# # Force relationships
	@constraint(model, [sc = 1:scenarios], buoyancy[:,sc] .== buoyancy_interpolated[sc])

	# Constraint (17) - Unsure what is correct
	@constraint(model, shear  .== cumulative_weight .- buoyancy)
	#@constraint(model, shear .== cumulative_weight .+ buoyancy)

	# @constraint(model, -100 <= sum(shear[i] * frame_length[i] for i in 1:n_positions-1) <= 100)

	@constraint(model,[sc = 1:scenarios], bending[1,sc] == 0)
	@constraint(model,[sc = 1:scenarios], bending[end,sc] == 0)
	# Constraint (19). 
	# Rasmus: Index for frame i i-1, because of the diff() function to calculate frame_length
	@constraint(model, [i = 2:n_positions, sc = 1:scenarios],
		bending[i,sc] == bending[i-1,sc] + (shear[i,sc] + shear[i-1,sc]) / 2 * frame_length[i-1]
	)
	stress_positions = [s.position for s in stress_limits]
	# Find the corresponding frame indices for these positions
	stress_frame_indices = [findmin(abs.(frame_positions .- pos))[2] for pos in stress_positions]

	# Stress limits. Constraint (18), (20)
	@constraint(model, [i in 1:length(stress_limits),sc = 1:scenarios],
		shear[stress_frame_indices[i],sc] + slack_shearMin[i,sc] >= shear_min[i])
	@constraint(model, [i in 1:length(stress_limits),sc = 1:scenarios],
		shear[stress_frame_indices[i],sc] <= shear_max[i] + slack_shearMax[i,sc])
	@constraint(model, [i in 1:length(stress_limits),sc = 1:scenarios],
		bending[stress_frame_indices[i],sc] <= bending_max[i] + slack_bendingMax[i,sc])

	# # Center of gravity constraints.
	@expression(model, lcg_ballast[sc = 1:scenarios],
		sum(vessel.ballast_tanks[t].lcg * ballast_volume[t,sc] for t in 1:n_ballast_tanks)
	)
	@expression(model, lcg_deadweight,
    sum(vessel.deadweight_tanks[t].lcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	# Left-hand side of constraint (11) & (12)
	@expression(model, lcg_total[sc = 1:scenarios],
    	lcg_cargo[sc] + lcg_ballast[sc] + lcg_deadweight + 
    vessel.lightship_lcg * vessel.lighship_weight
	)
	# Constraint (11)
	@constraint(model, lcg_max[sc = 1:scenarios],
		lcg_total[sc] <= sum(lcb[b] * displacement[b] * z_max[b,sc]
						 for b in draft_index_min:draft_index_max)
							+ slack_Lmax[sc]
	)
	# Constraint (12)
	@constraint(model, lcg_min[sc = 1:scenarios],
		lcg_total[sc] + slack_Lmin[sc] >= sum(lcb[b] * displacement[b] * z_min[b,sc]
						 for b in draft_index_min:draft_index_max)
	)
	# Left-hand side of constraint (13) & (14)
	@expression(model, tcg_ballast[sc = 1:scenarios],
		sum(vessel.ballast_tanks[t].tcg * ballast_volume[t,sc] for t in 1:n_ballast_tanks)
	)
	@expression(model, tcg_deadweight,
    sum(vessel.deadweight_tanks[t].tcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	@expression(model, tcg_total[sc = 1:scenarios],
		tcg_cargo[sc] + tcg_ballast[sc] + tcg_deadweight + 
		vessel.lightship_tcg * vessel.lighship_weight
	)
	# Constraint (13)
	@constraint(model, tcg_max[sc = 1:scenarios],
		tcg_total[sc] <= TCGmax * cumulative_weight[end] + slack_Tmax[sc]
	)
	# Constraint (14)
	@constraint(model, tcg_min[sc = 1:scenarios],
		tcg_total[sc] + slack_Tmin >= TCGmin * cumulative_weight[end]
	)
	# Some of the left-hand side of constraint (15) & (16)
	# Rasmus: vcg_slope relates to m_{hat}^{V,T} in the paper
	@expression(model, vcg_ballast[sc = 1:scenarios],
		sum(vcg_slope[t] * ballast_volume[t,sc] for t in 1:n_ballast_tanks)
	)
	@expression(model, vcg_deadweight,
    	sum(vessel.deadweight_tanks[t].vcg * vessel.deadweight_tanks[t].weight 
        for t in 1:length(vessel.deadweight_tanks))
	)
	# Changed 07/04
	@expression(model, vcg_total[sc = 1:scenarios],
		vcg_cargo[sc] + vcg_ballast[sc] + vcg_deadweight + 
		vessel.lightship_vcg * vessel.lighship_weight #vessel.lightship_vcg
	)

	# # Metacentric height constraint. Constraint (15) - NB DIFFERENT THAN IN ARTICLE
	@constraint(model, metacentric_height[sc = 1:scenarios],
		sum(gm[b] * displacement[b] * z_min[b,sc] for b in draft_index_min:draft_index_max) <=
		sum(metacenter[b] * displacement[b] * z_min[b,sc] for b in draft_index_min:draft_index_max) - vcg_total[sc]
		+ slack_Vmin[sc]
	)

	# KG constraint - NB DIFFERENT THAN IN ARTICLE
	@expression(model, kg_min_total[sc = 1:scenarios],
		sum(kg[b] * displacement[b] * z_max[b,sc] for b in draft_index_min:draft_index_max)
	)
	# Constraint (16)
	@constraint(model, kg_constraint[sc = 1:scenarios],
		kg_min_total[sc] + slack_Vmax[sc] >= vcg_total[sc]
	)
	return model
end
