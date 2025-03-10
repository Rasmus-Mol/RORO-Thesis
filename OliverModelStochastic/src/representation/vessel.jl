const DATA_DIR = joinpath(@__DIR__, "../../data/stability")

const REQUIRED_FILES = [
    "info.json",
    "ballast_tanks.csv",
    "deadweight.csv", 
    "weigt_dependencies.csv",
    "stress_limits.csv",
    "frame_positions.csv",
    "buoyancy.csv",
    "deck_weight_limits.csv"
]

@kwdef struct DeadWeight
    lcg::Float64
    tcg::Float64
    vcg::Float64
    weight::Float64
    tk_name::String # Rasmus: What is this?
    start_pos::Float64
    end_pos::Float64
    function DeadWeight(lcg, tcg, vcg, weight, tk_name, start_pos, end_pos)
        @assert start_pos < end_pos "Start position must be less than end position"
        @assert weight >= 0 "Weight must be non-negative"
        new(lcg, tcg, vcg, weight, tk_name, start_pos, end_pos)
    end
end

function Base.show(io::IO, t::DeadWeight)
    print(io, "DeadWeight(\"$(t.tk_name)\", pos=$(t.start_pos):$(t.end_pos), weight=$(round(t.weight, digits=1))t)")
end

@kwdef struct BallastTank
    max_vol::Float64
    start_pos::Float64
    end_pos::Float64
    lcg::Float64
    tcg::Float64
    min_vcg::Float64
    max_vcg::Float64
    tk_name::String
end

function Base.show(io::IO, t::BallastTank)
    print(io, "BallastTank(\"$(t.tk_name)\", pos=$(t.start_pos):$(t.end_pos), vol=$(round(t.max_vol, digits=1))m³)")
end

@kwdef struct WeightDependency
    displacement::Float64
    lcf::Float64  # Longitudinal Center of Floatation
    kmt::Float64  # Transverse Metacentric Height
    draft::Float64
    gm::Float64   # Metacentric Height
    kg::Float64   # Vertical Center of Gravity
    lcb::Float64  # Longitudinal Center of Buoyancy
end

function Base.show(io::IO, w::WeightDependency)
    print(io, "WeightDependency(disp=$(round(w.displacement, digits=1))t, draft=$(round(w.draft, digits=3))m)")
end

@kwdef struct StressLimit
    position::Float64
    bending_max::Float64
    bending_min::Union{Float64, Nothing}
    shear_max::Float64
    shear_min::Float64
end

function Base.show(io::IO, s::StressLimit)
    print(io, "StressLimit(pos=$(s.position)m, bend_max=$(round(s.bending_max))kNm)")
end


@kwdef struct FramePosition
    position::Float64
    lightship_weight_cumulative::Float64
end

function Base.show(io::IO, f::FramePosition)
    print(io, "Frame(pos=$(f.position)m, weight=$(round(f.lightship_weight_cumulative, digits=1))t)")
end

@kwdef struct BuoyancyPoint
    displacement::Float64
    cumulative_buoyancy::Float64
end

function Base.show(io::IO, b::BuoyancyPoint)
    print(io, "Buoyancy(disp=$(round(b.displacement))t, buoy=$(round(b.cumulative_buoyancy, digits=1))t)")
end

@kwdef struct Vessel
    # Basic vessel properties
    name::String
    IMO::Int
    lighship_weight::Float64
    lightship_lcg::Float64
    lightship_tcg::Float64
    lightship_vcg::Float64
    beam::Float64
    height_of_shaft_above_BL::Float64
    propeller_diameter::Float64
    
    # Tanks and loads
    ballast_tanks::Vector{BallastTank}
    deadweight_tanks::Vector{DeadWeight}

    decks::Vector{Deck}
    
    # Ship characteristics
    weight_dependencies::Vector{WeightDependency}
    stress_limits::Vector{StressLimit}
    frame_positions::Vector{FramePosition}
    buoyancy_points::Vector{BuoyancyPoint}
    buoyancy_displacement_weight_cumulative::Matrix{Float64}

    ballast_tank_frame_overlap::Matrix{Float64}  # [tank_idx, frame_idx]
    deadweight_tank_frame_overlap::Matrix{Float64}  # [tank_idx, frame_idx]
end

function Base.show(io::IO, v::Vessel)
    println(io, "Vessel \"$(v.name)\" (IMO: $(v.IMO))")
    println(io, "├─ Dimensions:")
    println(io, "│  ├─ Beam: $(v.beam)m")
    println(io, "│  └─ Propeller: ∅$(v.propeller_diameter)m at $(v.height_of_shaft_above_BL)m above BL")
    println(io, "├─ Lightship: $(round(v.lighship_weight, digits=1))t")
    println(io, "│  ├─ LCG: $(v.lightship_lcg)m")
    println(io, "│  ├─ TCG: $(v.lightship_tcg)m")
    println(io, "│  └─ VCG: $(v.lightship_vcg)m")
    println(io, "├─ Decks: $(length(v.decks)) decks")
    println(io, "├─ Ballast Tanks: $(length(v.ballast_tanks)) tanks")
    println(io, "├─ Deadweight Tanks: $(length(v.deadweight_tanks)) tanks")
    println(io, "├─ Weight Dependencies: $(length(v.weight_dependencies)) points")
    println(io, "├─ Stress Limits: $(length(v.stress_limits)) points")
    println(io, "├─ Frame Positions: $(length(v.frame_positions)) frames")
    print(io,   "└─ Buoyancy Points: $(length(v.buoyancy_points)) points")
end

# Rasums: bit unsure what this function does
function extract_buoyancy_displacement_weight_cumulative(displacement_cumulative_buoyancy, cumulative_buoyancy, displacement_draft)
    # Get unique displacement values and ensure non-empty
    x = unique(displacement_cumulative_buoyancy)
    N = length(x)
    if N == 0
        return zeros(length(displacement_draft), 1)
    end
    
    # Calculate points per displacement safely
    total_points = length(cumulative_buoyancy)
    if total_points == 0 || N == 0
        return zeros(length(displacement_draft), 1)
    end
    
    # Calculate points per displacement and ensure it's at least 1
    points_per_displacement = max(1, total_points ÷ N)
    
    # Create matrix to store values
    v = zeros(N, points_per_displacement)
    
    # Fill matrix with cumulative buoyancy values
    for n in 1:N
        disp_indices = findall(d -> isapprox(d, x[n], rtol=1e-10), displacement_cumulative_buoyancy)
        if !isempty(disp_indices)
            for m in 1:min(points_per_displacement, length(disp_indices))
                v[n, m] = cumulative_buoyancy[disp_indices[m]]
            end
        end
    end

    # Create interpolation functions for each position
    try
        li = [LinearInterpolation(x, v[:, m]) for m in 1:points_per_displacement]
        buoyancy_displacement_weight_cumulative = reduce(hcat, [li[m](displacement_draft) for m in eachindex(li)])
        return buoyancy_displacement_weight_cumulative
    catch e
        @warn "Failed to create interpolation, returning zeros" exception=e
        return zeros(length(displacement_draft), points_per_displacement)
    end
end

function create_deck_from_vessel_data(id::Integer, row::DataFrameRow)
    Deck(
        id = id,
        name = row.name,
        vcg = row.vcg,
        weight_limit = row.weight_limit,
        length = 0, #TODO Implement
        width = 0, #TODO Implement
        cargo_types = collect(values(CARGO_TYPES))  # Allow all cargo types by default
    )
end

# Rasmus: Calculates fraction of tank that overlaps with frame
function calculate_tank_frame_overlap(tank_start::Float64, tank_end::Float64, frame_start::Float64, frame_end::Float64)
    # If tank is completely outside frame
    if tank_end <= frame_start || tank_start >= frame_end
        return 0.0
    end
    
    # Calculate overlap
    overlap_start = max(tank_start, frame_start)
    overlap_end = min(tank_end, frame_end)
    overlap_length = overlap_end - overlap_start
    
    # Return percentage of tank that falls within this frame
    tank_length = tank_end - tank_start
    return overlap_length / tank_length
end
# Rasmus: Calculates overlap for all tanks and frames. Returns matrix of overlaps
function calculate_tank_overlaps(tanks::Vector{T}, frame_positions::Vector{FramePosition}) where T <: Union{BallastTank, DeadWeight}
    n_tanks = length(tanks)
    n_frames = length(frame_positions)
    overlaps = zeros(n_tanks, n_frames)
    
    for (t, tank) in enumerate(tanks)
        for p in 1:n_frames
            frame_start = p == 1 ? -Inf : frame_positions[p-1].position
            frame_end = frame_positions[p].position
            overlaps[t, p] = calculate_tank_frame_overlap(
                tank.start_pos,
                tank.end_pos,
                frame_start,
                frame_end
            )
        end
    end
    return overlaps
end


function Vessel(vessel_name::String)
    @info "Creating vessel object for '$vessel_name'"
    
    # Construct paths
    vessel_dir = joinpath(DATA_DIR, vessel_name)
    @debug "Looking for vessel data in: $vessel_dir"
    
    # Verify directory exists
    if !isdir(vessel_dir)
        @error "Vessel directory not found" vessel_dir
        throw(ArgumentError("Vessel data not found: $vessel_dir"))
    end
                     
    for file in REQUIRED_FILES
        path = joinpath(vessel_dir, file)
        if !isfile(path)
            @error "Required file missing" file path
            throw(ArgumentError("Required file not found: $path"))
        end
    end
    
    @info "All required files found, loading vessel data..."

    info_path = joinpath(vessel_dir, "info.json")
    ballast_tanks_path = joinpath(vessel_dir, "ballast_tanks.csv")
    deadweight_path = joinpath(vessel_dir, "deadweight.csv")
    weight_deps_path = joinpath(vessel_dir, "weigt_dependencies.csv")
    stress_limits_path = joinpath(vessel_dir, "stress_limits.csv")
    frame_positions_path = joinpath(vessel_dir, "frame_positions.csv")
    buoyancy_path = joinpath(vessel_dir, "buoyancy.csv")
    deck_weight_limits_path = joinpath(vessel_dir, "deck_weight_limits.csv")
    
    # Load vessel info
    vessel_data = open(info_path) do f
        JSON.parse(f)
    end
    
    # Load all CSV data
    tanks = CSV.read(ballast_tanks_path, DataFrame) |> 
           df -> [BallastTank(
               max_vol = row.max_vol,
               start_pos = row.start_pos,
               end_pos = row.end_pos,
               lcg = row.lcg,
               tcg = row.tcg,
               min_vcg = row.min_vcg,
               max_vcg = row.max_vcg,
               tk_name = row.tk_name
           ) for row in eachrow(df)]
    
    deadweights = CSV.read(deadweight_path, DataFrame) |>
                 df -> [DeadWeight(
                    lcg = row.lcg,
                    tcg = row.tcg,
                    vcg = row.vcg,
                    weight = row.weight,
                    tk_name = row.tk_name,
                    start_pos = row.start_pos,
                    end_pos = row.end_pos
                 ) for row in eachrow(df)]
    
    weight_deps = CSV.read(weight_deps_path, DataFrame) |>
                 df -> [WeightDependency(
                     displacement = row.displacement,
                     lcf = row.lcf,
                     kmt = row.kmt,
                     draft = row.draft,
                     gm = row.gm,
                     kg = row.kg,
                     lcb = row.lcb
                 ) for row in eachrow(df)]
    
    stress_limits = CSV.read(stress_limits_path, DataFrame) |>
                 df -> [StressLimit(
                     position = row.position,
                     bending_max = row.bending_max,
                     bending_min = ismissing(row.bending_min) ? nothing : row.bending_min,
                     shear_max = row.shear_max,
                     shear_min = row.shear_min
                 ) for row in eachrow(df)]
    
    frames = CSV.read(frame_positions_path, DataFrame) |>
            df -> [FramePosition(
                position = row.position,
                lightship_weight_cumulative = row.lightship_weight_cumulative
            ) for row in eachrow(df)]
    
    buoyancy = CSV.read(buoyancy_path, DataFrame) |>
               df -> [BuoyancyPoint(
                    displacement = row.displacement,
                    cumulative_buoyancy = row.cumulative_buoyancy
               ) for row in eachrow(df)]
    
    decks = CSV.read(deck_weight_limits_path, DataFrame) |>
            df -> [create_deck_from_vessel_data(i, row) 
                for (i, row) in enumerate(eachrow(df))]  



    displacement_cumulative_buoyancy = [b.displacement for b in buoyancy]
    cumulative_buoyancy = [b.cumulative_buoyancy for b in buoyancy]
    displacement_draft = [w.displacement for w in weight_deps]

    buoyancy_displacement_weight_cumulative = extract_buoyancy_displacement_weight_cumulative(displacement_cumulative_buoyancy, cumulative_buoyancy, displacement_draft)

    ballast_overlaps = calculate_tank_overlaps(tanks, frames)
    deadweight_overlaps = calculate_tank_overlaps(deadweights, frames)

    try
        name = String(vessel_data["name"])
        IMO = Int(vessel_data["IMO"])
        lighship_weight = Float64(vessel_data["lighship_weight"])
        lightship_lcg = Float64(vessel_data["lightship_lcg"])
        lightship_tcg = Float64(vessel_data["lightship_tcg"])
        lightship_vcg = Float64(vessel_data["lightship_vcg"])
        beam = Float64(vessel_data["beam"])
        height_of_shaft_above_BL = Float64(vessel_data["height_of_shaft_above_BL"])
        propeller_diameter = Float64(vessel_data["propeller_diameter"])

        # Use existing constructor with validated data
        vessel = Vessel(
            name = name,
            IMO = IMO,
            lighship_weight = lighship_weight,
            lightship_lcg = lightship_lcg,
            lightship_tcg = lightship_tcg,
            lightship_vcg = lightship_vcg,
            beam = beam,
            height_of_shaft_above_BL = height_of_shaft_above_BL,
            propeller_diameter = propeller_diameter,
            ballast_tanks = sort(tanks, by = x -> x.start_pos),
            deadweight_tanks = sort(deadweights, by = x -> x.start_pos),
            decks = sort(decks, by = x -> x.vcg),
            weight_dependencies = sort(weight_deps, by = x -> x.displacement),
            stress_limits = sort(stress_limits, by = x -> x.position),
            frame_positions = sort(frames, by = x -> x.position),
            buoyancy_points = sort(buoyancy, by = x -> x.displacement),
            buoyancy_displacement_weight_cumulative = buoyancy_displacement_weight_cumulative,
            ballast_tank_frame_overlap = ballast_overlaps,
            deadweight_tank_frame_overlap = deadweight_overlaps
        )
        @info "Successfully created vessel object" name=name IMO=IMO
        return vessel
    catch e
        @error "Failed to create vessel object" exception=e
        throw(ArgumentError("Error parsing vessel data: $e"))
    end
end

# Rasmus: Calculates stress at a given position.
# Rasmus: We only have data for stress at some point, so this function interpolates between points.
function interpolate_stress_limits(vessel::Vessel, position::Float64)
    stress_limits = vessel.stress_limits
    
    # Check if position is outside range
    pos_min = minimum(x -> x.position, stress_limits)
    pos_max = maximum(x -> x.position, stress_limits)
    
    if position < pos_min || position > pos_max
        throw(ArgumentError("Position $position is outside valid range [$pos_min, $pos_max]"))
    end
    
    # Find surrounding points
    idx_after = findfirst(x -> x.position >= position, stress_limits)
    idx_before = idx_after - 1
    
    # Handle exact match
    if stress_limits[idx_after].position == position
        limit = stress_limits[idx_after]
        # Rasmus: Don't understand why return statement is like this
        return (bending_max=limit.bending_max, 
                bending_min=limit.bending_min,
                shear_max=limit.shear_max,
                shear_min=limit.shear_min)
    end
    
    # Interpolate between points
    p1 = stress_limits[idx_before]
    p2 = stress_limits[idx_after]
    
    # Rasmus: t is the fraction of the distance between p1 and p2, from p1 to position
    t = (position - p1.position) / (p2.position - p1.position)
    
    bending_max = p1.bending_max + t * (p2.bending_max - p1.bending_max)
    bending_min = p1.bending_min + t * (p2.bending_min - p1.bending_min)
    shear_max = p1.shear_max + t * (p2.shear_max - p1.shear_max)
    shear_min = p1.shear_min + t * (p2.shear_min - p1.shear_min)
    # Rasmus: Don't understand why return statement is like this
    return (bending_max=bending_max, bending_min=bending_min, shear_max=shear_max, shear_min=shear_min)
end


function merge_frames(frame_positions::Vector{Float64}, frame_weights::Vector{Float64}, target_frames::Int)
    # Preserve first and last positions and weights
    first_pos = frame_positions[1]
    last_pos = frame_positions[end]
    first_weight = frame_weights[1]
    last_weight = frame_weights[end]
    
    # Remove first and last elements for merging
    interior_frames = frame_positions[2:end-1]
    interior_weights = frame_weights[2:end-1]
    adjusted_target = target_frames - 2
    
    # Calculate sections
    n_interior = length(interior_frames)
    quarter_idx = round(Int, n_interior * 0.25)
    three_quarter_idx = round(Int, n_interior * 0.75)
    
    # Calculate frame counts
    middle_frames = floor(Int, adjusted_target / 4)
    front_frames = floor(Int, (adjusted_target - middle_frames) / 2)
    back_frames = adjusted_target - front_frames - middle_frames
    
    # Split sections
    front_pos = interior_frames[1:quarter_idx]
    middle_pos = interior_frames[quarter_idx+1:three_quarter_idx]
    back_pos = interior_frames[three_quarter_idx+1:end]
    
    front_weights = interior_weights[1:quarter_idx]
    middle_weights = interior_weights[quarter_idx+1:three_quarter_idx]
    back_weights = interior_weights[three_quarter_idx+1:end]
    
    function merge_section(positions::Vector{Float64}, weights::Vector{Float64}, n_frames::Int)
        if length(positions) <= n_frames
            return positions, weights
        end
        
        indices = range(1, length(positions), length=n_frames)
        merged_pos = Float64[]
        merged_weights = Float64[]
        
        for i in 2:length(indices)
            start_idx = round(Int, indices[i-1])
            end_idx = round(Int, indices[i])
            push!(merged_pos, mean(positions[start_idx:end_idx]))
            push!(merged_weights, sum(weights[start_idx:end_idx]))
        end
        return merged_pos, merged_weights
    end
    
    # Merge sections
    merged_front_pos, merged_front_weights = merge_section(front_pos, front_weights, front_frames)
    merged_middle_pos, merged_middle_weights = merge_section(middle_pos, middle_weights, middle_frames)
    merged_back_pos, merged_back_weights = merge_section(back_pos, back_weights, back_frames)
    
    # Combine results
    merged_positions = vcat([first_pos], merged_front_pos, merged_middle_pos, merged_back_pos, [last_pos])
    merged_weights = vcat([first_weight], merged_front_weights, merged_middle_weights, merged_back_weights, [last_weight])
    
    # Calculate frame spacing
    frame_spacing = diff(merged_positions)
    
    return merged_positions, merged_weights, frame_spacing
end

# Rasmus: Generates n positions, with equal spacing.
function merge_positions(positions::Vector{Float64}, n::Int)
    if n >= length(positions)
        return positions  # Return original if requested size is larger
    end
    
    # Get total length of ship
    total_length = positions[end] - positions[1]
    
    # Calculate new section length
    section_length = total_length / (n-1)
    
    # Create new positions array
    new_positions = zeros(n)
    new_positions[1] = positions[1]  # Keep first position
    new_positions[end] = positions[end]  # Keep last position
    
    # Fill in intermediate positions with equal spacing
    for i in 2:(n-1)
        new_positions[i] = positions[1] + (i-1) * section_length
    end
    
    return new_positions
end

# Rasmus: Generates new frame positions and buoyancy points with equal spacing, length = target_points
# does the same for buoyancy points
function merge_buoyancy_points_and_frames(vessel::Vessel, target_points::Int)
    if target_points >= length(vessel.frame_positions)
        return vessel.frame_positions, vessel.buoyancy_points
    end
    
    # Get frame positions and weights
    frame_pos = [f.position for f in vessel.frame_positions]
    frame_weights = [f.lightship_weight_cumulative for f in vessel.frame_positions]
    
    # Calculate new positions with equal spacing
    new_positions = range(frame_pos[1], frame_pos[end], length=target_points)
    
    # Interpolate weights at new positions
    weight_interp = LinearInterpolation(frame_pos, frame_weights)
    new_weights = weight_interp.(new_positions)
    
    # Create new frame positions
    new_frame_positions = [FramePosition(
        position = pos,
        lightship_weight_cumulative = weight
    ) for (pos, weight) in zip(new_positions, new_weights)]
    
    # Handle buoyancy points for each displacement
    unique_displacements = unique([b.displacement for b in vessel.buoyancy_points])
    new_buoyancy_points = Vector{BuoyancyPoint}()
    
    for displacement in unique_displacements
        # Get points for this displacement
        disp_points = filter(b -> b.displacement == displacement, vessel.buoyancy_points)
        positions = [vessel.frame_positions[i].position for i in 1:length(vessel.frame_positions)]
        buoyancy_values = [p.cumulative_buoyancy for p in disp_points]
        
        # Create interpolation for buoyancy values
        buoyancy_interp = LinearInterpolation(positions, buoyancy_values)
        
        # Interpolate at new positions
        for pos in new_positions
            push!(new_buoyancy_points, BuoyancyPoint(
                displacement = displacement,
                cumulative_buoyancy = buoyancy_interp(pos)
            ))
        end
    end
    
    return new_frame_positions, new_buoyancy_points
end

# Rasmus: Simplifies vessel by reducing number of frames and buoyancy points
# so the hydrostatic table is smaller
function simplify_vessel(vessel::Vessel; target_points::Int=20)
    @info "Simplifying vessel with target points: $target_points"
    
    # Merge frame positions and buoyancy points together
    new_frame_positions, new_buoyancy_points = merge_buoyancy_points_and_frames(vessel, target_points)
    
    # Create new matrix for buoyancy displacement weight cumulative
    displacement_draft = [w.displacement for w in vessel.weight_dependencies][1:2:end] # <- Rasmus: What is iteration?
    displacement_cumulative_buoyancy = [b.displacement for b in new_buoyancy_points]
    cumulative_buoyancy = [b.cumulative_buoyancy for b in new_buoyancy_points]
    
    new_buoyancy_matrix = extract_buoyancy_displacement_weight_cumulative(
        displacement_cumulative_buoyancy,
        cumulative_buoyancy,
        displacement_draft
    )

    new_ballast_overlaps = calculate_tank_overlaps(vessel.ballast_tanks, new_frame_positions)
    new_deadweight_overlaps = calculate_tank_overlaps(vessel.deadweight_tanks, new_frame_positions)
    
    # Create new vessel with simplified data
    return Vessel(
        name=vessel.name,
        IMO=vessel.IMO,
        lighship_weight=vessel.lighship_weight,
        lightship_lcg=vessel.lightship_lcg,
        lightship_tcg=vessel.lightship_tcg,
        lightship_vcg=vessel.lightship_vcg,
        beam=vessel.beam,
        height_of_shaft_above_BL=vessel.height_of_shaft_above_BL,
        propeller_diameter=vessel.propeller_diameter,
        ballast_tanks=vessel.ballast_tanks,
        deadweight_tanks=vessel.deadweight_tanks,
        decks=vessel.decks,
        weight_dependencies=[w for w in vessel.weight_dependencies][1:2:end],
        stress_limits=vessel.stress_limits,
        frame_positions=new_frame_positions,
        buoyancy_points=new_buoyancy_points,
        buoyancy_displacement_weight_cumulative=new_buoyancy_matrix,
        ballast_tank_frame_overlap = new_ballast_overlaps,
        deadweight_tank_frame_overlap = new_deadweight_overlaps
    )
end

# vessel = Vessel("finlandia")
# println("Original vessel:")
# println("- Frame positions: $(length(vessel.frame_positions))")
# println("- Buoyancy points: $(length(vessel.buoyancy_points))")

# simplified_vessel = simplify_vessel(vessel, target_points=25)
# println("\nSimplified vessel:")
# println("- Frame positions: $(length(simplified_vessel.frame_positions))")
# println("- Buoyancy points: $(length(simplified_vessel.buoyancy_points))")

# # Verify relationship is maintained
# n_frames = length(simplified_vessel.frame_positions)
# n_unique_displacements = length(unique([b.displacement for b in simplified_vessel.buoyancy_points]))
# points_per_displacement = length(simplified_vessel.buoyancy_points) ÷ n_unique_displacements

# println("\nVerification:")
# println("- Points per displacement: $points_per_displacement")
# println("- Number of frames: $n_frames")
# @assert points_per_displacement == n_frames "Relationship between frames and buoyancy points maintained"