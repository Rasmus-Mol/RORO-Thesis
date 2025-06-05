# Rasmus: Creates an array with slope vcg for all tanks
function calculate_vcg_slopes(vessel::Vessel)
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

# Rasmus: Add stability constraints to the model
function add_stability!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)
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
	@expression(model, vcg_ballast,
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


# Add stability constraints with slack variables to the model
# Slack added at: Centers of gravity
function add_stability_slack!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo, 
		slack_Vmax,slack_Vmin,slack_Tmin,slack_Tmax,slack_Lmin,slack_Lmax, # slack for center of gravity
		slack_shear1,slack_shear2,slack_shearMin,slack_shearMax,slack_bendingMax,slack_ballast_tanks) # slack for stress and bending

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
		ballast_volume[t] <= vessel.ballast_tanks[t].max_vol * ρ + slack_ballast_tanks[t]
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
	# FIXED
	@expression(model, buoyancy_interpolated,
		sum(vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_min[b] + 
		vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_max[b] 
			for b in draft_index_min:draft_index_max)./2
	)
	# # Force relationships
	@constraint(model, buoyancy .== buoyancy_interpolated)

	# Constraint (17) - First is correct
	#@constraint(model, shear .+ slack_shear2 .== cumulative_weight .- buoyancy .+ slack_shear1)
	@constraint(model, shear .== cumulative_weight .- buoyancy)

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
		shear[stress_frame_indices[i]] + slack_shearMin[i] >= shear_min[i])
	@constraint(model, [i in 1:length(stress_limits)],
		shear[stress_frame_indices[i]] <= shear_max[i] + slack_shearMax[i])
	@constraint(model, [i in 1:length(stress_limits)],
		bending[stress_frame_indices[i]] <= bending_max[i] + slack_bendingMax[i])

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
							+ slack_Lmax
	)
	# Constraint (12)
	@constraint(model, lcg_min,
		lcg_total + slack_Lmin >= sum(lcb[b] * displacement[b] * z_min[b]
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
		tcg_total <= TCGmax * cumulative_weight[end] + slack_Tmax
	)
	# Constraint (14)
	@constraint(model, tcg_min,
		tcg_total + slack_Tmin >= TCGmin * cumulative_weight[end]
	)
	# Some of the left-hand side of constraint (15) & (16)
	# Rasmus: vcg_slope is m_{hat}^{V,T} in paper
	@expression(model, vcg_ballast,
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
		+ slack_Vmin
	)

	# KG constraint
	@expression(model, kg_min_total,
		sum(kg[b] * displacement[b] * z_max[b] for b in draft_index_min:draft_index_max)
	)
	# Constraint (16)
	@constraint(model, kg_constraint,
		kg_min_total + slack_Vmax >= vcg_total
	)

	return model
end

# Allows for no ballast water
function add_stability_no_ballast!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo)
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

	@variable(model, ballast_volume[1:n_ballast_tanks] >= 0)  # Ballast volumes
	@variable(model, cumulative_weight[1:n_positions] >= 0)
	@variable(model, z_min[draft_index_min:draft_index_max], Bin)
	@variable(model, z_max[draft_index_min:draft_index_max], Bin)

	# Force ballast tanks to be empty
	@constraint(model, sum(ballast_volume) == 0
	)

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
	@expression(model, vcg_ballast,
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


# Add stability constraints with slack variables to the model
# No water allowed in tanks
function add_stability_slack_no_water!(vessel::Vessel, model, pos_weight_cargo, lcg_cargo, tcg_cargo, vcg_cargo,
    slack_Vmax, slack_Vmin, slack_Tmin, slack_Tmax, slack_Lmin, slack_Lmax, # slack for center of gravity
    slack_shear1, slack_shear2, slack_shearMin, slack_shearMax, slack_bendingMax, slack_ballast_tanks) # slack for stress and bending

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

    @variable(model, ballast_volume[1:n_ballast_tanks] >= 0)  # Ballast volumes
    @variable(model, cumulative_weight[1:n_positions] >= 0)
    @variable(model, z_min[draft_index_min:draft_index_max], Bin)
    @variable(model, z_max[draft_index_min:draft_index_max], Bin)

    # Force variables
    @variable(model, buoyancy[1:n_positions])
    @variable(model, shear[1:n_positions])
    @variable(model, bending[1:n_positions])

    # New
    @constraint(model, [t in 1:n_ballast_tanks], ballast_volume[t] == 0)

    # Ballast water constraints. Constraint (9)
    @constraint(
        model,
        [t = 1:n_ballast_tanks],
        ballast_volume[t] <= vessel.ballast_tanks[t].max_vol * ρ + slack_ballast_tanks[t]
    )

    # Draft index constraints. Constraint (3), (4), (5)
    @constraint(model, sum(z_min) == 1)
    @constraint(model, sum(z_max) == 1)
    @constraint(model, [b = draft_index_min:draft_index_max-1], z_min[b] - z_max[b+1] == 0)

    # Calculate ballast tank weights per position. Constraint (10)
    @expression(model, pos_weight_tank[p=1:n_positions],
        sum(ballast_volume[t] * vessel.ballast_tank_frame_overlap[t, p]
            for t in 1:n_ballast_tanks)
    )

    @expression(model, pos_weight_deadweight[p=1:n_positions],
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

	# New slack variables
	@variable(model, slack_dis_min >= 0)
	@variable(model, slack_dis_max >= 0)
    # Draft constraints. Constraint (6)
    @constraint(model, dis_min,
        sum(displacement[b] * z_min[b] for b in draft_index_min:draft_index_max) -
        cumulative_weight[end] <= 0 + slack_dis_min
    )
    # Constraint (7)
    @constraint(model, dis_max,
        sum(displacement[b] * z_max[b] for b in draft_index_min:draft_index_max) -
        cumulative_weight[end] + slack_dis_max >= 0
    )

    # Buoyancy calculations using vessel's buoyancy matrix
    # Rasmus: I think this should be constraint (8) but it doesn't look right
    # FIXED
    @expression(model, buoyancy_interpolated,
        sum(vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_min[b] +
            vessel.buoyancy_displacement_weight_cumulative[b, :] .* z_max[b]
            for b in draft_index_min:draft_index_max) ./ 2
    )
    # # Force relationships
    @constraint(model, buoyancy .== buoyancy_interpolated)

    # Constraint (17) - First is correct
    #@constraint(model, shear .+ slack_shear2 .== cumulative_weight .- buoyancy .+ slack_shear1)
    @constraint(model, shear .== cumulative_weight .- buoyancy)

    # @constraint(model, -100 <= sum(shear[i] * frame_length[i] for i in 1:n_positions-1) <= 100)

    @constraint(model, bending[1] == 0)
    @constraint(model, bending[end] == 0)
	# new
	@variable(model, slack_bending_pos[1:n_positions] >= 0)
	@variable(model, slack_bending_neg[1:n_positions] >= 0)
    # Constraint (19). 
    # Rasmus: Index for frame i i-1, because of the diff() function to calculate frame_length
    @constraint(model, [i = 2:n_positions],
        bending[i] == bending[i-1] + (shear[i] + shear[i-1]) / 2 * frame_length[i-1] + slack_bending_pos[i] - slack_bending_neg[i]
    )

    stress_positions = [s.position for s in stress_limits]
    # Find the corresponding frame indices for these positions
    stress_frame_indices = [findmin(abs.(frame_positions .- pos))[2] for pos in stress_positions]

    # Stress limits. Constraint (18), (20)
    @constraint(model, [i in 1:length(stress_limits)],
        shear[stress_frame_indices[i]] + slack_shearMin[i] >= shear_min[i])
    @constraint(model, [i in 1:length(stress_limits)],
        shear[stress_frame_indices[i]] <= shear_max[i] + slack_shearMax[i])
    @constraint(model, [i in 1:length(stress_limits)],
        bending[stress_frame_indices[i]] <= bending_max[i] + slack_bendingMax[i])

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
                     +
                     slack_Lmax
    )
    # Constraint (12)
    @constraint(model, lcg_min,
        lcg_total + slack_Lmin >= sum(lcb[b] * displacement[b] * z_min[b]
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
        tcg_total <= TCGmax * cumulative_weight[end] + slack_Tmax
    )
    # Constraint (14)
    @constraint(model, tcg_min,
        tcg_total + slack_Tmin >= TCGmin * cumulative_weight[end]
    )
    # Some of the left-hand side of constraint (15) & (16)
    # Rasmus: vcg_slope is m_{hat}^{V,T} in paper
    @expression(model, vcg_ballast,
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
        +
        slack_Vmin
    )

    # KG constraint
    @expression(model, kg_min_total,
        sum(kg[b] * displacement[b] * z_max[b] for b in draft_index_min:draft_index_max)
    )
    # Constraint (16)
    @constraint(model, kg_constraint,
        kg_min_total + slack_Vmax >= vcg_total
    )

    return model
end
