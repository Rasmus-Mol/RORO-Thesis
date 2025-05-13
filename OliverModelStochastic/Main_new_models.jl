# Script to run on HPC
push!(LOAD_PATH, pwd())
# Doens't work without this, dont know why
include("src/StowagePlannerStochastic.jl")
using .StowagePlannerStochastic
using JuMP

parse_index = parse(Int, ARGS[1]) # Jobindex input


problem_det = load_data("finlandia", "no_cars_medium_100_haz_eq_0.1", "hazardous")

# Folder name for results - date and hour
HPC_folder = "Finlandia_"*Dates.format(now(), "dd_mm_HH")*"New_Models"
 # Describe tests if necessary
extra_info = "Change this"

# First job index - create problem 
if parse_index == 1
    # Creates Deterministic problem and model
    model_det = create_model(problem_det)
    set_silent(model_det) # removes terminal output
    set_time_limit_sec(model_det, 60 * 60) # 1 hour
    optimize!(model_det)
    solution_det = extract_solution(problem_det, model_det)
    # Save Solution for deterministic problem
    write_solution(solution_det,"Finlandia_deterministic","Deterministic_Solution",HPC_folder)
end

# Check if folder and file has been created, otherwise create
file_check = "Results/"*HPC_folder*"HPC_data.json"
if !isfile(file_check)
    # Save scenario and number of unknown weights
    write_HPC_data(repetitions, scenarios, n_cargo_unknownweight, time_limit, HPC_folder, extra_info)
end

number_of_cargo = [80,90,100,length(problem_det.cargo)]
n = number_of_cargo[parse_index] # number of cargo for current job
# Run tests
carogc = CawrgoCollection([shuffle!(problem.cargo.items)[i] for i in 1:n])
ntypes = [length(filter(x -> x.cargo_type_id == i, problem.cargoc.items)) for i in 1:4]
r_cargo = random_cargocollection(ntypes)
for i in 1:repetitions
    # From # original problem
    r_plan1, not_s_r = random_stowage_plan(cargoc, problem_det.slots)
    mo = create_random_stowageplan_model(r_plan1, not_s_r, cargoc, vessel, slots, false)
    set_silent(mo) # removes terminal output
    set_time_limit_sec(mo, 60)
    optimize!(mo)
    pro = StowageProblem(vessel = vessel,slots = slots,cargo = cargoc,name = det_problem.name)
    sol = extract_solution(pro,mo)random_cargocollection(ntypes::Vector{Int}, shuffle::Bool = true)
    # From random cargocollection
    r_plan2, not_s_r2 = random_stowage_plan(r_cargo, problem_det.slots)
    random_cargocollection(ntypes::Vector{Int}, shuffle::Bool = true)
end
create_model_stochastic_cargo_fraction(problem::StochasticStowageProblem, fraction::Float64)

random_cargocollection(ntypes::Vector{Int}, shuffle::Bool = true)



# choose n first elements and shuffle them
n = 117
cargoc = CargoCollection(shuffle!([cargo.items[i] for i in 1:n]))
splan, not_s = random_stowage_plan(cargoc, slots)
println("Number of cargoes loaded: ", n-length(not_s),"/",n)
#splan2, not_s2 = random_stowage_plan(sort_cargocollection(cargoc), slots)
# New random cargocollection
ntypes = [length(filter(x -> x.cargo_type_id == i, cargoc.items)) for i in 1:4]
rcargoc = random_cargocollection(ntypes)
splan_r, not_s_r = random_stowage_plan(rcargoc, slots)
println("Number of cargoes loaded: ", n-length(not_s_r),"/",n)
#splan_r2, not_s_r2 = random_stowage_plan(sort_cargocollection(rcargoc), slots)
