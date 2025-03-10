@kwdef struct LMCargo
    wght::Float64
    h::Float64
    l::Float64
    b::Float64
    lcg::Float64
    tcg::Float64
    dVCG::Float64
    name::String
end

@kwdef struct LoadMaster
    apiroot::String
    apiversion::Union{Integer, AbstractFloat}
    headers::Dict
end
LoadMaster(apiroot; apiversion=2.0, headers=Dict()) = LoadMaster(apiroot, apiversion, headers)

@kwdef struct Item
    status::String
    val::String
    unit::String
    code::String
end

@kwdef struct TankSimple
    TankIndex:: Int64
    volume:: Float64
end

@kwdef struct TankFull
    fsm::Float64
    shortTankName::String
    colorHex::String
    numberOfTankGroup::Int
    prcVolume::Float64
    productName::String
    tcg::Float64
    vcg::Float64
    lcg::Float64
    fullNameTankGroup::String
    volume::Float64
    weight::Float64
    density::Float64
    tankIndex::Int
    maxVolume::Float64
    shortNameTankGroup::String
    fullTankName::String
end

@kwdef struct LMVessel
    name::String
    length::Float64
    order::Int64
end

function LMVessel(exp::Dict{String,Any})
    name = exp["name"]
    length = parse(Float64, replace(exp["length"], "," =>"."))
    order = parse(Int64, exp["order"])
    LMVessel(name=name, length=length, order=order)
end
