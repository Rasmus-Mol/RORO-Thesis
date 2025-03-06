function lmget(lm::LoadMaster, endpoint::String; kwargs...)
    apiuri = uri(lm, endpoint, kwargs)
    apiheaders = headers(lm, Dict("Content-Type" => "application/json"))

    try
        response = HTTP.get(apiuri, apiheaders)
        return JSON.parse(String(response.body))
    catch e
        throw(e)
    end
end

function lmpost(lm::LoadMaster, endpoint::String, body::String)
    apiuri = uri(lm, endpoint)
    apiheaders = headers(lm, Dict("Content-Type" => "application/json"))

    try
        response = HTTP.post(apiuri, apiheaders, body)
        return JSON.parse(String(response.body))
    catch e
        throw(e)
    end
end

function lmdelete(lm, endpoint; kwargs...)
    apiuri = uri(lm, endpoint, kwargs)
    apiheaders = headers(lm, Dict("Content-Type" => "application/json"))

    try
        response = HTTP.delete(apiuri, apiheaders)
        return JSON.parse(String(response.body))
    catch e
        throw(e)
    end
end

function deleteallcargo(lm::LoadMaster)
    @info "Deleting all cargo"
    try
        endpoint = "drycargos/ALL"
        result = lmdelete(lm, endpoint)
        return result
    catch e
        if isa(e, HTTP.ExceptionRequest.StatusError) && e.status == 404
            return missing
        end
        throw(e)
    end
end

function loadcargo(lm::LoadMaster, cargos::Vector{LMCargo})
    @info "Loading $(length(cargos)) cargo(s)"
    try
        endpoint = "drycargos"
        cargos_json = JSON.json(cargos)
        result = lmpost(lm, endpoint, cargos_json)
        return result
    catch e
        throw(e)
    end
end

function gettanks(lm::LoadMaster)
    @info "Getting tanks"
    try
        endpoint = "tanks"
        result = lmget(lm, endpoint)
        return result
    catch e
        throw(e)
    end
end

function loadtanks(lm::LoadMaster, tanks::Vector{TankSimple})
    @info "Loading tanks*"
    try
        endpoint = "tanks"
        ballast_json = JSON.json(tanks)
        result = lmpost(lm, endpoint, ballast_json)
        return result
    catch e
        throw(e)
    end
end

function resettanks(lm::LoadMaster, n_tanks::Int64)
    tanks = Vector{TankSimple}(undef, n_tanks)
    for i in 1:n_tanks
        tanks[i] = TankSimple(i, 0)
    end

    @info "Resetting tanks"
    loadtanks(lm, tanks)
end

function parse_tanks(tanks::Dict{String, Any})
    tanks_data = tanks["tanks"]
    tanks = [TankFull(Float64(t["fsm"]), 
                      t["shortTankName"], 
                      t["colorHex"],
                      t["numberOfTankGroup"], 
                       t["prcVolume"],
                      t["productName"], 
                      t["tcg"], 
                      t["vcg"], 
                      t["lcg"], 
                      t["fullNameTankGroup"], 
                      t["volume"], 
                      t["weight"], 
                      t["density"],
                      t["tankIndex"], 
                      t["maxVolume"],
                      t["shortNameTankGroup"], 
                      t["fullTankName"]) for t in tanks_data]
    return tanks
end

function getresults(lm::LoadMaster)
    @info "Getting results"
    try
        endpoint = "results"
        result = lmget(lm, endpoint)
        return result
    catch e
        throw(e)
    end
end

function parse_result(result::Dict{String, Any})
    items = [Item(item["status"], item["val"], item["unit"], item["code"]) for item in result["itemList"]]
    return items
end

function uri(lm::LoadMaster, endpoint="", query=missing)
    u = URI("$(lm.apiroot)/$(endpoint)")
    !ismissing(query) && return URI(u; query=query)
    u
end

headers(lm::LoadMaster, custom_headers::AbstractDict) = merge(lm.headers, custom_headers)

function struct_to_json(type::Type{T}) where T
    return JSON.json(type)
end

function getcurrentvessel(lm::LoadMaster)
    @info "Getting current vessel"
    try
        endpoint = "vessel"
        result = lmget(lm, endpoint)
        return LMVessel(result)
    catch e
        if isa(e, HTTP.ExceptionRequest.StatusError) && e.status == 404
            return missing
        end
        throw(e)
    end
end

function resetvessel(lm::LoadMaster)
    @info "Resetting vessel"
    deleteallcargo(lm)
    resettanks(lm, 21)
end

function upload_solution(lm::LoadMaster, sol)
    # Convert cargo placements to LMCargo format
    cargos = [LMCargo(
                  wght = c.weight,
                  h = 2.0,  # Assuming 2D placement
                  l = c.length,
                  b = c.width,
                  lcg = c.lcg,
                  tcg = c.tcg,
                  dVCG = c.vcg,
                  name = "Cargo_$(c.id)"
              ) for c in sol.cargo]

    # Convert ballast volumes to TankSimple objects 
    tanks = [TankSimple(i, vol) for (i, vol) in enumerate(sol.ballast_volume)]

    # Upload to API
    @info "Uploading solution with $(length(cargos)) cargos and $(length(tanks)) tanks"
    cargo_result = loadcargo(lm, cargos)
    tanks_result = loadtanks(lm, tanks)

    return (cargo_result = cargo_result, tanks_result = tanks_result)
end

# lm = LoadMaster("http://localhost:5000")
# v = getcurrentvessel(lm)
# resetvessel(lm)

# upload_solution(lm, solution)

# result = getresults(lm)
# parse_result(result)