# Add metadata struct
@kwdef struct InstanceConfig
    name::String
    capacity_ratio::Float64
    refrigerated::Bool
    hazardous::Bool
    weight_profile::String
end

# Update metadata struct to include config
@kwdef struct InstanceMetadata
    name::String
    description::String
    total_cargo::Int
    total_weight::Float64
    cargo_counts::Dict{String,Int}
    weight_profile::String
    capacity_ratio::Float64
    hazardous_enabled::Bool
    refrigerated_enabled::Bool
    generation_timestamp::DateTime
    config::InstanceConfig
end

function load_cargo_from_loadmaster(filepath::String)
    df = DataFrame(XLSX.readtable(filepath, "data"))

    cargos = Vector{Cargo}()
    for row in eachrow(df)
        cargo = Cargo(
            id = string(row.G_RefNo),
            cargo_type = SpecialCargo,  # Default to special cargo for loadmaster data
            weight = row.G_Weight,
            loading_port = 1,  # Default values since loadmaster doesn't specify ports
            discharge_port = 2,
            priority = 1,
            requires_lashing = false,
            requires_ventilation = false,
            hazardous = 0,
            refers = false
        )
        push!(cargos, cargo)
    end

    return CargoCollection(cargos)
end

function load_instance(filepath::String)
    # Read and parse JSON
    instance = JSON3.read(read(filepath, String))
    
    # Parse metadata
    cfg = instance.config
    config = InstanceConfig(
        name=cfg.name,
        capacity_ratio=cfg.capacity_ratio,
        refrigerated=cfg.refrigerated,
        hazardous=cfg.hazardous,
        weight_profile=cfg.weight_profile
    )
    
    # Parse metadata
    meta = instance.metadata
    metadata = InstanceMetadata(
        name=meta.name,
        description=meta.description,
        total_cargo=meta.total_cargo,
        total_weight=meta.total_weight,
        cargo_counts=Dict(String(k) => v for (k,v) in meta.cargo_counts),
        weight_profile=meta.weight_profile,
        capacity_ratio=meta.capacity_ratio,
        hazardous_enabled=meta.hazardous_enabled,
        refrigerated_enabled=meta.refrigerated_enabled,
        generation_timestamp=DateTime(meta.generation_timestamp),
        config=config
    )
    
    # Convert cargo items
    cargos = Vector{Cargo}()
    for (i, item) in enumerate(instance.cargo)
        cargo_type = if item.type == "Trailer"
            Truck
        elseif item.type == "Secu Box"
            Secu
        elseif item.type == "Machine"
            Machine
        else
            Car  # Default fallback
        end

        cargo = Cargo(
            id = item.id,
            cargo_type_id = Int(cargo_type),
            weight = item.weight,
            loading_port = 1,  # Default values since not specified in JSON
            discharge_port = 2,
            priority = 1,
            requires_lashing = false,
            requires_ventilation = false,
            hazardous = item.hazard_class,
            refers = item.refrigerated
        )
        push!(cargos, cargo)
    end
    
    return CargoCollection(cargos), metadata
end

