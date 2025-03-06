module StowagePlanner

using Base: @kwdef
import Base: length, getindex, lastindex, keys, eachindex, iterate, eltype, firstindex

using Random
using StructArrays
using JSON3, JSON
using DataFrames
using CSV, XLSX
using Dates
using Interpolations
using UnPack
using JuMP, HiGHS, GLPK
using Plots
using HTTP
using URIs

include("representation/cargo.jl")
include("representation/deck.jl")
include("representation/slot.jl")
include("representation/vessel.jl")
include("representation/instance.jl")
include("representation/problem.jl")

include("model/base_model.jl")
include("model/hazardous.jl")
include("model/stability.jl")

include("utils/helpers.jl")
include("solution.jl")
include("plots/solution.jl")

include("loadmaster/types.jl")
include("loadmaster/api.jl")


# Makes these functions available outside the module
export load_data, create_model, solve_model, extract_solution, plot_solution
export StowageProblem

export LoadMaster, getcurrentvessel, resetvessel, upload_solution, getresults, parse_result, uri, headers, struct_to_json

end # module StowagePlanner
