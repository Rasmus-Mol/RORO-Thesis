@kwdef struct CargoPlacement
    id::Int
    cargo_type_id::Int
    deck::Int
    lcg::Float64
    tcg::Float64
    vcg::Float64
    weight::Float64
    length::Float64
    width::Float64
    height::Float64
    haz_class::Int
end

@kwdef struct Solution

    gap::Float64
    status::Any
    objective::Float64
    time::Float64

    cargo_weight::Float64
    total_weight::Float64
    ballast_weight::Float64

    area_utilization::Float64 # calculation is not correct
    cargo::Vector{CargoPlacement}
    n_cargo_total::Int = 0
    n_cargo_loaded::Int = 0

    shear_force::Vector{Float64}
    bending_moment::Vector{Float64}
    ballast_volume::Vector{Float64}
    lcg::Float64
    tcg::Float64
    vcg::Float64

    n_variables::Int = 0
    n_binary_variables::Int = 0
    n_constraints::Int = 0
    model_size::Int = 0  # Memory size in bytes
    solver_name::String = ""
    solver_iterations::Int = 0
    solver_nodes::Int = 0

end

function Base.show(io::IO, sol::Solution)
    println(io, "Solution Status:")
    println(io, "├─ Gap:              $(round(sol.gap * 100, digits=3))%")
    println(io, "├─ Status:           $(sol.status)")
    println(io, "├─ Objective:        $(round(sol.objective, digits=3))")
    println(io, "├─ Solve Time:       $(round(sol.time, digits=3)) seconds")
    println(io, "├─ Model Statistics:")
    println(io, "│  ├─ Variables:     $(sol.n_variables) ($(sol.n_binary_variables) binary)")
    println(io, "│  ├─ Constraints:   $(sol.n_constraints)")
    println(io, "│  └─ Size:          $(round(sol.model_size/1024/1024, digits=2)) MB")
    println(io, "├─ Solver Details:")
    println(io, "│  ├─ Name:          $(sol.solver_name)")
    println(io, "│  ├─ Iterations:    $(sol.solver_iterations)")
    println(io, "│  └─ Nodes:         $(sol.solver_nodes)")
    println(io, "├─ Cargo Loading:")
    println(io, "│  ├─ Loaded:        $(sol.n_cargo_loaded)/$(sol.n_cargo_total)")
    println(io, "│  ├─ Weight:        $(round(sol.cargo_weight, digits=3)) t")
    println(io, "│  └─ Utilization:   $(round(sol.area_utilization * 100, digits=3))%")
    println(io, "├─ Weights:")
    println(io, "│  ├─ Cargo:         $(round(sol.cargo_weight, digits=3)) t")
    println(io, "│  ├─ Ballast:       $(round(sum(sol.ballast_volume), digits=3)) t")
    println(io, "│  └─ Total:         $(round(sol.total_weight, digits=3)) t")
    println(io, "└─ Center of Gravity:")
    println(io, "   ├─ LCG:           $(round(sol.lcg, digits=3)) m")
    println(io, "   ├─ TCG:           $(round(sol.tcg, digits=3)) m")
    println(io, "   └─ VCG:           $(round(sol.vcg, digits=3)) m")
end

function extract_solution(problem, model)
    # First check if solution exists
    status = termination_status(model)
    if status ∈ [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
        return Solution(
            gap = Inf,
            status = status,
            objective = Inf,
            time = solve_time(model),
            cargo_weight = 0.0,
            total_weight = 0.0,
            ballast_weight = 0.0, 
            area_utilization = 0.0,
            cargo = CargoPlacement[],
            n_cargo_total = length(problem.cargo),
            n_cargo_loaded = 0,
            shear_force = Float64[],
            bending_moment = Float64[],
            ballast_volume = Float64[],
            lcg = 0.0,
            tcg = 0.0, 
            vcg = 0.0,
            n_variables = num_variables(model),
            n_binary_variables = count(is_binary, all_variables(model)),
            n_constraints = num_constraints(model; count_variable_in_set_constraints=true),
            model_size = Base.summarysize(model),
            solver_name = string(JuMP.solver_name(model)),
            solver_iterations = 0,
            solver_nodes = 0
        )
    end

    # Extract solution values
    cs_val = value.(model[:cs])
    ballast_val = reverse(value.(model[:ballast_volume]))
    weight_val = value.(model[:weight])
    
    # Calculate statistics
    cargo_weight = sum(weight_val)
    ballast_weight = sum(ballast_val)
    total_weight = value.(model[:cumulative_weight][end])

    # Get placed cargo
    placements = CargoPlacement[] # Rasmus: Same as Vector{CargoPlacement}()
    n_cargo_loaded = 0
    
    for c in 1:length(problem.cargo), s in 1:length(problem.slots)
        if cs_val[c,s] > 0.5  # Binary variable threshold
            n_cargo_loaded += 1
            slot = problem.slots[s]
            cargo = problem.cargo[c]
            push!(placements, CargoPlacement(
                id = cargo.id,
                cargo_type_id = cargo.cargo_type_id,
                deck = slot.deck_id,
                lcg = slot.lcg,
                tcg = slot.tcg,
                vcg = slot.vcg,
                weight = cargo.weight,
                length = get_length(cargo),
                width = get_width(cargo),
                height = get_height(cargo),
                haz_class = cargo.hazardous
            ))
        end
    end

    # Calculate center of gravity
    lcg = value(model[:lcg_total]) / total_weight
    tcg = value(model[:tcg_total]) / total_weight
    vcg = value(model[:vcg_total]) / total_weight

    # Calculate area utilization
    total_cargo_area = sum(p.length * p.width for p in placements)
    total_vessel_area = sum(
        slot.length * slot.width
        for slot in problem.slots
    ) # Rasmus: Don't the slots have overlapping areas, so this is larger than the actual area?
    area_utilization = total_cargo_area / total_vessel_area

    return Solution(
        gap = relative_gap(model),
        status = status,
        objective = objective_value(model),
        time = solve_time(model),
        cargo_weight = cargo_weight,
        total_weight = total_weight,
        ballast_weight = ballast_weight,
        area_utilization = area_utilization,
        cargo = placements,
        n_cargo_total = length(problem.cargo),
        n_cargo_loaded = n_cargo_loaded,
        shear_force = [], # These are not in the base model
        bending_moment = [], # These are not in the base model
        ballast_volume = collect(ballast_val),
        lcg = lcg,
        tcg = tcg,
        vcg = vcg,
        n_variables = num_variables(model),
        n_binary_variables = count(is_binary, all_variables(model)),
        n_constraints = num_constraints(model; count_variable_in_set_constraints=true),
        model_size = Base.summarysize(model),
        solver_name = string(JuMP.solver_name(model)),
        solver_iterations = 0,
        solver_nodes = 0
    )
end