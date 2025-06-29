# Write data to JSON files

# Save problem
# Save 1 cargo collection
function write_problem(problem::StowageProblem,foldername::String,filename::String, HPC_folder::String)
    # Creates HPC_folder for results for this HPC instance
    if !isdir("Results/"*HPC_folder)
        mkdir("Results/"*HPC_folder)
    end
    # Creates folder for results
    if !isdir("Results/"*HPC_folder*"/"*foldername)
        mkdir("Results/"*HPC_folder*"/"*foldername)
    end
    cargo = [i for i in problem.cargo]
    open(joinpath("Results",HPC_folder,foldername,filename*".json"), "w") do file
        JSON.print(file, cargo, 4)  # Pretty-print with indentation
    end
    info = [problem.name,problem.timestamp]
    open(joinpath("Results",HPC_folder,foldername,filename*"_info"*".json"), "w") do file
        JSON.print(file, info, 4)  # Pretty-print with indentation
    end
end
# Saves the cargo collection of a stochastic problem
function write_problem_stochastic(problem::StochasticStowageProblem,foldername::String,filename::String,HPC_folder::String)
    # Creates folder for results
    if !isdir("Results/"*HPC_folder*"/"*foldername)
        mkdir("Results/"*HPC_folder*"/"*foldername)
    end
    for i in 1:problem.scenarios
        cargo = [j for j in problem.cargo.items[i]]
        open(joinpath("Results",HPC_folder,foldername,filename*"_scenario_"*string(i)*".json"), "w") do file
            JSON.print(file, cargo, 4)  # Pretty-print with indentation
        end
    end
    info = [problem.name,problem.unknown_weights,problem.known_weights,problem.scenarios,problem.probability,problem.timestamp]
    open(joinpath("Results",HPC_folder,foldername,filename*"_info"*".json"), "w") do file
        JSON.print(file, info, 4)  
    end
end
# Save solutions
# Saves deterministic solution
function write_solution(solution::Solution,foldername::String,filename::String,HPC_folder::String)
    # Creates HPC_folder for results for this HPC instance
    if !isdir("Results/"*HPC_folder)
        mkdir("Results/"*HPC_folder)
    end
    if !isdir("Results/"*HPC_folder*"/"*foldername)
        mkdir("Results/"*HPC_folder*"/"*foldername)
    end
    placements = [cargo for cargo in solution.cargo]
    info = [solution.gap,solution.status, solution.objective, 
        solution.time, solution.cargo_weight, solution.total_weight,
        solution.ballast_weight, solution.area_utilization,solution.cs,
        solution.n_cargo_total, solution.n_cargo_loaded, solution.shear_force,
        solution.bending_moment, solution.ballast_volume, solution.lcg,
        solution.tcg, solution.vcg,solution.n_variables, solution.n_binary_variables,
        solution.n_constraints, solution.model_size, solution.solver_name,
        solution.solver_iterations, solution.solver_nodes]
    open(joinpath("Results",HPC_folder,foldername,filename*"_placements"*".json"), "w") do file
        JSON.print(file, placements, 4)  # Pretty-print with indentation
    end
    open(joinpath("Results",HPC_folder,foldername,filename*"_info"*".json"), "w") do file
        JSON.print(file, info, 4)  # Pretty-print with indentation
    end
end
function write_solution_stochastic(solution::SolutionStochastic,foldername::String,filename::String,HPC_folder::String)
    # Creates HPC_folder for results for this HPC instance
    if !isdir("Results/"*HPC_folder)
        mkdir("Results/"*HPC_folder)
    end
    if !isdir("Results/"*HPC_folder*"/"*foldername)
        mkdir("Results/"*HPC_folder*"/"*foldername)
    end
    placements = [cargo for cargo in solution.cargo]
    forces = [force for force in solution.forces]
    info = [solution.gap,solution.status, solution.objective,
            solution.objective_scenario, 
        solution.time, solution.area_utilization,solution.cs,
        solution.n_cargo_total, solution.n_cargo_loaded, 
        solution.n_variables, solution.n_binary_variables,
        solution.n_constraints, solution.model_size, solution.solver_name,
        solution.solver_iterations, solution.solver_nodes]
    open(joinpath("Results",HPC_folder,foldername,filename*"_placements"*".json"), "w") do file
        JSON.print(file, placements, 4)  # Pretty-print with indentation
    end
    open(joinpath("Results",HPC_folder,foldername,filename*"_forces"*".json"), "w") do file
        JSON.print(file, forces, 4)
    end
    open(joinpath("Results",HPC_folder,foldername,filename*"_info"*".json"), "w") do file
        JSON.print(file, info, 4)
    end
end

# Get data back from JSON files
# Returns deterministic problem
function get_deterministic_problem(foldername::String,filename::String,HPC_folder::String,problemname1,problemname2,problemname3)
    open(joinpath("Results",HPC_folder,foldername,filename*".json"), "r") do file
        cargo = JSON3.read(read(file, String), Vector{Cargo})
        # load ship data
        problem = load_data(problemname1, problemname2, problemname3)
        open(joinpath("Results",HPC_folder,foldername,filename*"_info"*".json"), "r") do file
            info = JSON3.read(read(file, String), Vector{Any})
            name = info[1]
            timestamp = info[2]
            det_pro = StowageProblem(
                vessel = problem.vessel,
                slots = problem.slots,
                cargo = CargoCollection(cargo),
                name = name,
                timestamp = DateTime(timestamp)
            )
            return det_pro
        end
    end
end
# Returns stochastic problem
function get_stochastic_problem(foldername::String,filename::String,HPC_folder::String,problemname1,problemname2,problemname3)
    open(joinpath("Results",HPC_folder,foldername,filename*"_info"*".json"), "r") do file
        info = JSON3.read(read(file, String), Vector{Any})
        name = info[1]
        unknown_weights = info[2]
        known_weights = info[3]
        scenarios = info[4]
        probability = info[5]
        timestamp = info[6]
        # load ship data
        problem = load_data(problemname1, problemname2, problemname3)

        cargo = Vector{CargoCollection}()
        for i in 1:scenarios
            open(joinpath("Results",HPC_folder,foldername,filename*"_scenario_"*string(i)*".json"), "r") do file
                cargo_c = CargoCollection(JSON3.read(read(file, String), Vector{Cargo}))
                push!(cargo, cargo_c)
            end
        end
        stochastic_problem = StochasticStowageProblem(
            vessel = problem.vessel,
            slots = problem.slots,
            cargo = CargoCollectionScenarios(cargo),
            unknown_weights = unknown_weights,
            known_weights = known_weights,
            scenarios = scenarios,
            probability = probability,
            name = name,
            timestamp = DateTime(timestamp)
        )
        return stochastic_problem
    end
end

# Get deterministic solution
function get_solution_deterministic(foldername::String,filename::String,HPC_folder::String)
    open(joinpath("Results",HPC_folder,foldername,filename*"_info"*".json"), "r") do file
        info = JSON3.read(read(file, String), Vector{Any})
        gap = info[1]
        if isnothing(gap)
            gap = Inf # no solution
        end
        status = info[2]
        objective = info[3]
        if isnothing(objective)
            objective = Inf # no solution
        end
        time = info[4]
        cargo_weight = info[5]
        total_weight = info[6]
        ballast_weight = info[7]
        area_utilization = info[8]
        cs = info[9]
        # convert to correct type
        if length(cs)>0
            cs_M = Matrix{Float64}(undef, length(cs[1]), length(cs))
            for i in 1:length(cs)
                for j in 1:length(cs[1])
                    cs_M[j,i] = cs[i][j]
                end
            end
        else # no solution
            cs_M = Matrix{Float64}(undef, 0, 0) # empty matrix
            gap = Inf
            objective = Inf
        end
        n_cargo_total = info[10]
        n_cargo_loaded = info[11]
        shear_force = info[12]
        bending_moment = info[13]
        ballast_volume = info[14]
        lcg = info[15]
        tcg = info[16]
        vcg = info[17]
        n_variables = info[18]
        n_binary_variables = info[19]
        n_constraints = info[20]
        model_size = info[21]
        solver_name = info[22]
        solver_iterations = info[23]
        solver_nodes = info[24]
        placements = open(joinpath("Results",HPC_folder,foldername,filename*"_placements"*".json"), "r") do file
            JSON3.read(read(file, String), Vector{CargoPlacement})
        end
        solution = Solution(
            gap = gap,
            status = status,
            objective = objective,
            time = time,
            cargo_weight = cargo_weight,
            total_weight = total_weight,
            ballast_weight = ballast_weight,
            cargo = placements,
            area_utilization = area_utilization,
            cs = cs_M,
            n_cargo_total = n_cargo_total,
            n_cargo_loaded = n_cargo_loaded,
            shear_force = shear_force,
            bending_moment = bending_moment,
            ballast_volume = ballast_volume,
            lcg = lcg,
            tcg = tcg,
            vcg = vcg,
            n_variables = n_variables,
            n_binary_variables = n_binary_variables,
            n_constraints = n_constraints,
            model_size = model_size,
            solver_name = solver_name,
            solver_iterations = solver_iterations,
            solver_nodes = solver_nodes
        )
        return solution
    end
end
# Get Stochastic solution
function get_solution_stochastic(foldername::String,filename::String,HPC_folder::String)
    open(joinpath("Results",HPC_folder,foldername,filename*"_info"*".json"), "r") do file
        info = JSON3.read(read(file, String), Vector{Any})
        gap = info[1]
        if isnothing(gap)
            gap = Inf # no solution
        end
        status = info[2]
        objective = info[3]
        if isnothing(objective)
            objective = Inf # no solution
        end
        objective_scenario = info[4]
        time = info[5]
        area_utilization = info[6]
        cs = info[7]
        # convert to correct type
        cs_M = Matrix{Float64}(undef, length(cs[1]), length(cs))
        for i in 1:length(cs)
            for j in 1:length(cs[1])
                cs_M[j,i] = cs[i][j]
            end
        end
        n_cargo_total = info[8]
        n_cargo_loaded = info[9]
        n_variables = info[10]
        n_binary_variables = info[11]
        n_constraints = info[12]
        model_size = info[13]
        solver_name = info[14]
        solver_iterations = info[15]
        solver_nodes = info[16]
        placements = open(joinpath("Results",HPC_folder,foldername,filename*"_placements"*".json"), "r") do file
            JSON3.read(read(file, String), Vector{CargoPlacementStochastic})
        end
        forces = open(joinpath("Results",HPC_folder,foldername,filename*"_forces"*".json"), "r") do file
            JSON3.read(read(file, String), Vector{WeightAndForceStochastic})
        end

        solution = SolutionStochastic(
            gap = gap,
            status = status,
            objective = objective,
            objective_scenario = objective_scenario,
            time = time,
            n_cargo_total = n_cargo_total,
            n_cargo_loaded = n_cargo_loaded,
            cargo = placements,
            area_utilization = area_utilization,
            cs = cs_M,
            forces = forces,
            n_variables = n_variables,
            n_binary_variables = n_binary_variables,
            n_constraints = n_constraints,
            model_size = model_size,
            solver_name = solver_name,
            solver_iterations = solver_iterations,
            solver_nodes = solver_nodes
        )
        return solution
    end
end

# Save HPC input and create folder for the following results
function write_HPC_data(repetitions::Int64, scenarios::Vector, n_unknown::Vector,time_limit,HPC_folder::String, info::String = nothing)
    data = [repetitions, scenarios, n_unknown, time_limit, info]
    # Creates folder for results
    if !isdir("Results/"*HPC_folder)
        mkdir("Results/"*HPC_folder)
    end
    open(joinpath("Results",HPC_folder,"HPC_data.json"), "w") do file
        JSON.print(file, data, 4)  # Pretty-print with indentation
    end
end
# Load HPC input
function get_HPC_data(HPC_folder::String)
    open(joinpath("Results",HPC_folder,"HPC_data.json"), "r") do file
        info = JSON3.read(read(file, String),Vector{Any})
        repetitions = info[1]
        scenarios = info[2]
        n_unknown = info[3]
        time_limit = info[4]
        if length(info) > 4
            info_string = info[5]
        else
            info_string = "No extra info"
        end
        return repetitions, scenarios, n_unknown, time_limit, info_string
    end
end

# Save slack variables
function write_slack(HPC_folder::String, foldername::String, filename::String, model)
    slack = [value.(model[:slack_deck]), value.(model[:slack_Vmax]),value.(model[:slack_Vmin]),
            value.(model[:slack_Tmin]), value.(model[:slack_Tmax]), 
            value.(model[:slack_Lmin]), value.(model[:slack_Lmax]),
            #value.(model[:slack_shear1]), value.(model[:slack_shear2]), 
            value.(model[:slack_shearMin]),
            value.(model[:slack_shearMax]), value.(model[:slack_bendingMax]),
             value.(model[:slack_ballast_tanks])
            ]   
    # Creates HPC_folder for results for this HPC instance
    if !isdir("Results/"*HPC_folder)
        mkdir("Results/"*HPC_folder)
    end
    if !isdir("Results/"*HPC_folder*"/"*foldername)
        mkdir("Results/"*HPC_folder*"/"*foldername)
    end
    open(joinpath("Results",HPC_folder,foldername,filename*"_slack.json"), "w") do file
        JSON.print(file, slack, 4)
    end
end
# Load slack variables
function get_slack(foldername::String,filename::String,HPC_folder::String)
    open(joinpath("Results",HPC_folder,foldername,filename*"_slack.json"), "r") do file
        slack = JSON3.read(read(file, String),Vector{Any})
        slack_deck = slack[1]
        slack_Vmax = slack[2]
        slack_Vmin = slack[3]
        slack_Tmin = slack[4]
        slack_Tmax = slack[5]
        slack_Lmin = slack[6]
        slack_Lmax = slack[7]
        #slack_shear1 = slack[8]
        #slack_shear2 = slack[9]
        slack_shearMin = slack[8]
        slack_shearMax = slack[9]
        slack_bending = slack[10]
        slack_ballast_tanks = slack[11]
        slack_var = [slack_deck, slack_Vmin, slack_Vmax, slack_Tmin, slack_Tmax,
                slack_Lmin, slack_Lmax, #slack_shear1, slack_shear2,
                slack_shearMin, slack_shearMax, slack_bending, slack_ballast_tanks]
        return slack_var
        #slack_deck, slack_Vmin, slack_Vmax, slack_Tmin, slack_Tmax,
         #       slack_Lmin, slack_Lmax, #slack_shear1, slack_shear2,
         #       slack_shearMin, slack_shearMax, slack_bending, slack_ballast_tanks
    end
end

# To save stability of solution tests
function save_stability_of_solution(HPC_folder::String, foldername::String, filename::String, feasible)
    if !isdir("Results/"*HPC_folder)
        mkdir("Results/"*HPC_folder)
    end
    if !isdir("Results/"*HPC_folder*"/"*foldername)
        mkdir("Results/"*HPC_folder*"/"*foldername)
    end
    m = [collect(feasible[i, :]) for i in 1:size(feasible, 1)]
    open(joinpath("Results",HPC_folder,foldername,filename*".json"), "w") do file
        JSON.print(file, m, 4)
    end
end 

function get_stability_of_solution(foldername::String, filename::String, HPC_folder::String)
    open(joinpath("Results",HPC_folder,foldername,filename*".json"), "r") do file
        feasible = JSON3.read(read(file, String), Vector{Vector{Float64}})
        return feasible
    end
end
