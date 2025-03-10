@kwdef struct Deck
    id::Integer
    name::String
    vcg::Float64    # Vertical Center of Gravity (height from baseline)
    weight_limit::Float64  # maximum weight capacity in tons
    length::Float64  # meters
    width::Float64   # meters
    cargo_types::Vector{CargoTypeInfo} = CargoTypeInfo[]  # allowed cargo types
end

# Helper functions - fix the field names to match the struct
get_vcg(deck::Deck) = deck.vcg
get_weight_limit(deck::Deck) = deck.weight_limit
get_dimensions(deck::Deck) = (deck.length, deck.width)

# Pretty printing for single deck
function Base.show(io::IO, ::MIME"text/plain", deck::Deck)
    println(io, "Deck(ID: $(deck.id))")
    println(io, "├─ Name: $(deck.name)")
    println(io, "├─ VCG: $(deck.vcg)m")
    println(io, "├─ Dimensions: $(deck.length)m × $(deck.width)m")
    println(io, "├─ Weight Limit: $(deck.weight_limit) tons")
    println(io, "└─ Allowed Cargo Types:")
    for (i, cargo_type) in enumerate(deck.cargo_types)
        prefix = i == length(deck.cargo_types) ? "   └─ " : "   ├─ "
        println(io, prefix, cargo_type.name)
    end
end

function Base.show(io::IO, deck::Deck)
    print(io, "Deck(name=$(deck.name), vcg=$(deck.vcg), weight_limit=$(deck.weight_limit)t)")
end

struct DeckCollection
    items::StructArray{Deck}
    total_capacity::Float64
    
    function DeckCollection(items::Vector{Deck})
        struct_items = StructArray(items)
        # Rasmus: Unsure about max_load here, should it be weight_limit?
        total_capacity = sum(d -> d.max_load, items)
        new(struct_items, total_capacity)
    end
end
# Rasmus: These are methods
length(collection::DeckCollection) = length(collection.items)
getindex(collection::DeckCollection, i::Int64) = collection.items[i]

# Filtering functions
function get_decks_by_cargo_type(collection::DeckCollection, cargo_type::CargoTypeInfo)
    filter(deck -> cargo_type in deck.cargo_types, collection.items)
end

function get_decks_by_height(collection::DeckCollection, min_height::Float64)
    filter(deck -> deck.height >= min_height, collection.items)
end

# Pretty printing for collection
function Base.show(io::IO, ::MIME"text/plain", collection::DeckCollection)
    println(io, "DeckCollection:")
    println(io, "├─ Total decks: $(length(collection.items))")
    println(io, "├─ Total capacity: $(collection.total_capacity) tons")
    println(io, "└─ Decks:")
    for (i, deck) in enumerate(collection.items)
        prefix = i == length(collection.items) ? "   └─ " : "   ├─ "
        println(io, prefix, deck)
    end
end

# Test data generation
function generate_test_deck(id::Integer)
    # Rasmus: Returns array with values from CARGO_TYPES
    cargo_types = collect(values(CARGO_TYPES))
    # Rasmus: Randomly select cargo types that are allowed
    allowed_types = rand(cargo_types, rand(1:length(cargo_types)))
    
    Deck(
        id=id,
        name="Deck $(id)",
        cargo_types=allowed_types,
        height=id * 3.0,  # 3m between decks
        max_load=rand(1000.0:100.0:5000.0),
        length=rand(100.0:10.0:300.0),
        width=rand(20.0:2.0:40.0)
    )
end

# function generate_test_collection(n::Int)
#     decks = [generate_test_deck(i) for i in 1:n]
#     DeckCollection(decks)
# end
