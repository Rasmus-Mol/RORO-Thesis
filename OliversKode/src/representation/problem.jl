@kwdef struct StowageProblem
    # Core components
    vessel::Vessel
    slots::SlotCollection
    cargo::CargoCollection    
    # Problem metadata
    name::String
    timestamp::DateTime = now()
end

function load_data(vessel_name::String, instance_name::String, instance_type::String)
    # Load vessel data
    vessel = Vessel(vessel_name)
    vessel = simplify_vessel(vessel, target_points=50)

    
    # Load slots from data directory
    slots = load_all_slots(joinpath(@__DIR__, "../../data/slots", vessel_name))
    
    # Load cargo and metadata from default instance
    instance_path = joinpath(@__DIR__, "../../instances/", instance_type, vessel_name, "$instance_name.json")
    cargo, metadata = load_instance(instance_path)
    
    return StowageProblem(
        vessel = vessel,
        slots = slots,
        cargo = cargo,
        name = metadata.name
    )
end