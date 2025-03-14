
@kwdef struct CargoTypeInfo
    id::Int
    name::String
    length::Float64  # meters
    width::Float64   # meters
    height::Float64  # meters
end

const CARGO_TYPES = Dict(
    "Car" => CargoTypeInfo(id=2, name="Cars", length=4.0, width=1.5, height=1.6),
    "Truck" => CargoTypeInfo(id=1, name="Truck", length=13.6, width=2.5, height=4.0),
    "Machine" => CargoTypeInfo(id=3, name="Machine", length=4.6, width=4.5, height=3.0),
    "Secu" => CargoTypeInfo(id=4, name="Secu", length=13.8, width=3.6, height=4.0),
)

@enum CargoType Unknown Truck Car Machine Secu

@kwdef struct Cargo
    id::Int
    cargo_type_id::Int64
    weight::Float64  # tons
    loading_port::Int64
    discharge_port::Int64
    priority::Int64 = 1
    requires_lashing::Bool = false
    requires_ventilation::Bool = false
    hazardous::Int = 0
    refers::Bool = false
end

# Update the get_type_info function to work with cargo_type_id instead of CargoType
function get_type_info(cargo::Cargo)
    for (name, info) in CARGO_TYPES
        if info.id == cargo.cargo_type_id
            return info
        end
    end
    error("No cargo type info found for cargo_type_id: $(cargo.cargo_type_id)")
end

# Rasmus: These are methods
# Update the dimension getter functions
get_length(cargo::Cargo) = get_type_info(cargo).length
get_width(cargo::Cargo) = get_type_info(cargo).width
get_height(cargo::Cargo) = get_type_info(cargo).height

# Update the show function to display the cargo type name
function Base.show(io::IO, ::MIME"text/plain", cargo::Cargo)
    type_info = get_type_info(cargo)
    println(io, "Cargo(ID: $(cargo.id))")
    println(io, "├─ Type: $(type_info.name)")
    println(io, "├─ Dimensions: $(get_length(cargo))m × $(get_width(cargo))m × $(get_height(cargo))m")
    println(io, "├─ Weight: $(cargo.weight) tons")
    println(io, "├─ Route: Port $(cargo.loading_port) → Port $(cargo.discharge_port)")
    println(io, "├─ Priority Level: $(cargo.priority)")
    println(io, "├─ Requires Lashing: $(cargo.requires_lashing)")
    println(io, "├─ Requires Ventilation: $(cargo.requires_ventilation)")
    println(io, "├─ Hazardous Level: $(cargo.hazardous > 0 ? "$(cargo.hazardous) ☣️" : "None")")  # conditional hazardous level
    println(io, "└─ Requires Refrigeration: $(cargo.refers ? "Yes ❄️" : "No")")  # display refrigeration requirement with emoji
end

function Base.show(io::IO, cargo::Cargo)
    print(io, "Cargo($(cargo.id), $(cargo.cargo_type_id), $(get_length(cargo))×$(get_width(cargo))×$(get_height(cargo))m, Hazardous: $(cargo.hazardous > 0 ? ", $(cargo.hazardous) ☣️, " : "None"), Refers: $(cargo.refers ? "Yes ❄️" : "No"))")
end

struct CargoCollection
    items::StructArray{Cargo}
    total_weight::Float64
    cargo_types::Array{Int, 1}
    
    # Inner constructor to calculate total weight and unique cargo types
    function CargoCollection(items::Vector{Cargo})
        struct_items = StructArray(items)
        total_weight = sum(struct_items.weight)
        cargo_types = unique(item.cargo_type_id for item in items)
        new(struct_items, total_weight, cargo_types)
    end
end

# A struct for each CargoCollection for each scenario
struct CargoCollectionScenarios
    items::StructArray{CargoCollection}
    total_weights::Vector{Float64}
    # constructor
    function CargoCollectionScenarios(items::Vector{CargoCollection})
        struct_items = StructArray(items)
        total_weights = [item.total_weight for item in items]
        new(struct_items, total_weights)
    end
end

length(collection::CargoCollection) = length(collection.items)

getindex(collection::CargoCollection, i::Int64) = collection.items[i]

# Define iteration interface
iterate(collection::CargoCollection) = iterate(collection.items)
iterate(collection::CargoCollection, state) = iterate(collection.items, state)

# Optional: Also implement eltype and firstindex, lastindex for completeness
eltype(::Type{CargoCollection}) = Cargo
firstindex(collection::CargoCollection) = firstindex(collection.items)
lastindex(collection::CargoCollection) = lastindex(collection.items)

keys(collection::CargoCollection) = 1:length(collection)
eachindex(collection::CargoCollection) = keys(collection)

function Base.filter(f::Function, collection::CargoCollection)
    filtered_items = filter(f, collection.items)
    CargoCollection(collect(filtered_items))
end

height(collection::CargoCollection, i::Int64) = get_height(collection[i])
Base.getproperty(collection::CargoCollection, i::Int64, ::Val{:height}) = height(collection, i)

# Get width of cargo at index i
width(collection::CargoCollection, i::Int64) = get_width(collection[i])
Base.getproperty(collection::CargoCollection, i::Int64, ::Val{:width}) = width(collection, i)

# Get length of cargo at index i
length(collection::CargoCollection, i::Int64) = get_length(collection[i])
Base.getproperty(collection::CargoCollection, i::Int64, ::Val{:length}) = length(collection, i)

# Update get_cargo_by_type function
function get_cargo_by_type(collection::CargoCollection, type::CargoType)
    type_id = CARGO_TYPES[string(type)].id
    filter(cargo -> cargo.cargo_type_id == type_id, collection.items)
end

function get_cargo_by_port(collection::CargoCollection, port::Int64)
    filter(cargo -> cargo.loading_port == port || cargo.discharge_port == port, collection.items)
end

# Update cargo type distribution function
function get_cargo_type_distribution(collection::CargoCollection)
    type_counts = Dict{Int, Int}()
    for cargo in collection.items
        type_counts[cargo.cargo_type_id] = get(type_counts, cargo.cargo_type_id, 0) + 1
    end
    type_counts
end

# Pretty printing
function Base.show(io::IO, ::MIME"text/plain", collection::CargoCollection)
    println(io, "CargoCollection:")
    println(io, "├─ Total items: $(length(collection.items))")
    println(io, "├─ Total weight: $(collection.total_weight) tons")
    println(io, "├─ Cargo types:")
    
    type_dist = get_cargo_type_distribution(collection)
    total = length(collection.items)
    for (i, (type, count)) in enumerate(sort(collect(type_dist)))
        percentage = round(count/total * 100, digits=1)
        prefix = i == length(type_dist) ? "│  └─ " : "│  ├─ "
        println(io, prefix, "$type: $count ($percentage%)")
    end
    
    println(io, "└─ Items:")
    for (i, cargo) in enumerate(collection.items)
        prefix = i == length(collection.items) ? "   └─ " : "   ├─ "
        println(io, prefix, cargo)
    end
end

#function generate_random_cargo(id::String)
function generate_random_cargo(id::Int)
    cargo_types = [Car, Truck, Machine, Secu]
    cargo_type = rand(cargo_types)
    type_info = CARGO_TYPES[string(cargo_type)]
    
    # Weight based on cargo type and volume
    volume = type_info.length * type_info.width * type_info.height
    
    weight = if cargo_type == Car
        rand(1.5:0.1:3.0)
    elseif cargo_type == Truck
        rand(5.0:0.5:40.0)
    else
        volume * rand(0.2:0.01:0.5) # density factor
    end
    
    ports = 1:5
    loading_port = rand(ports)
    discharge_port = rand(setdiff(ports, loading_port))
    
    Cargo(
        id=id,
        cargo_type_id=type_info.id,  # Use the id from type_info
        weight=weight,
        loading_port=loading_port,
        discharge_port=discharge_port,
        priority=rand(1:3),
        requires_lashing=rand(Bool),
        requires_ventilation=rand(Bool),
        hazardous=rand(0:3),
        refers=rand(Bool)
    )
end

function generate_test_collection(n::Int)
    #cargos = [generate_random_cargo("CARGO_$(lpad(i,3,'0'))") for i in 1:n]
    cargos = [generate_random_cargo(i) for i in 1:n]
    CargoCollection(cargos)
end
