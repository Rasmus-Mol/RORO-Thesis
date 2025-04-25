@kwdef struct Slot
    id::Int  # Unique identifier
    loadmaster_id::Int  # Original ID from loadmaster
    deck_id::Int
    cargo_type_id::Int
    lcg::Float64
    tcg::Float64
    vcg::Float64
    refrigerated::Bool
    length::Float64
    width::Float64
    on_deck::Bool     
end

# Helper functions to get slot properties
get_length(slot::Slot) = slot.length
get_width(slot::Slot) = slot.width

# Pretty printing for single slot
function Base.show(io::IO, ::MIME"text/plain", slot::Slot)
    println(io, "Slot(ID: $(slot.id), LoadmasterID: $(slot.loadmaster_id))")
    println(io, "├─ Deck: $(slot.deck_id)")
    println(io, "├─ Cargo Type: $(slot.cargo_type_id)")
    println(io, "├─ Dimensions: $(slot.length)m × $(slot.width)m")
    println(io, "├─ Position: LCG=$(slot.lcg), TCG=$(slot.tcg), VCG=$(slot.vcg)")
    println(io, "├─ Refrigerated: $(slot.refrigerated ? "Yes ❄️" : "No")")
    println(io, "└─ Location: $(slot.on_deck ? "On Deck 🔲" : "Below Deck 🔳")")
end
# Rasmus: Pretty sure it overwrites the previous show function. 
# IO is the stream, i.e. terminal or file
function Base.show(io::IO, slot::Slot)
    print(io, "Slot($(slot.id)[LM:$(slot.loadmaster_id)], Deck $(slot.deck_id), $(slot.length)×$(slot.width)m, $(slot.on_deck ? "On" : "Below") Deck)")
end

function do_slots_overlap(slot1::Slot, slot2::Slot, margin::Float64=0.0)
    # Only check overlap if slots are on the same deck
    if slot1.deck_id != slot2.deck_id
        return false
    end

    # Calculate slot boundaries using LCG (center) and dimensions, including margin
    slot1_left = slot1.lcg - slot1.length/2 - margin
    slot1_right = slot1.lcg + slot1.length/2 + margin
    slot1_front = slot1.tcg - slot1.width/2 - margin
    slot1_back = slot1.tcg + slot1.width/2 + margin

    slot2_left = slot2.lcg - slot2.length/2 - margin
    slot2_right = slot2.lcg + slot2.length/2 + margin
    slot2_front = slot2.tcg - slot2.width/2 - margin
    slot2_back = slot2.tcg + slot2.width/2 + margin

    # Check if rectangles overlap
    horizontal_overlap = !(slot1_right < slot2_left || slot2_right < slot1_left)
    vertical_overlap = !(slot1_back < slot2_front || slot2_back < slot1_front)

    return horizontal_overlap && vertical_overlap
end

function find_overlapping_slots(slots::AbstractArray{Slot, 1}, margin::Float64=0.0)
    n = length(slots)
    overlap_matrix = falses(n, n)
    
    for i in 1:n
        for j in (i+1):n
            if do_slots_overlap(slots[i], slots[j], margin)
                overlap_matrix[i,j] = overlap_matrix[j,i] = true
            end
        end
    end
    
    return overlap_matrix
end

struct SlotCollection
    items::StructArray{Slot}
    total_slots::Int
    overlap_matrix::Array{Bool, 2}
    
    function SlotCollection(items::Vector{Slot})
        struct_items = StructArray(items)
        total_slots = length(items)
        new(struct_items, total_slots, find_overlapping_slots(items, 0.1))
    end
end

length(collection::SlotCollection) = length(collection.items)
getindex(collection::SlotCollection, i::Int64) = collection.items[i]
firstindex(collection::SlotCollection) = firstindex(collection.items)
lastindex(collection::SlotCollection) = lastindex(collection.items)
keys(collection::SlotCollection) = Base.OneTo(length(collection))
eachindex(collection::SlotCollection) = keys(collection)
Base.iterate(collection::SlotCollection) = iterate(collection.items)
Base.iterate(collection::SlotCollection, state) = iterate(collection.items, state)

# Rasmus: Removes items for which function f returns false
function Base.filter(f::Function, collection::SlotCollection)
    filtered_items = filter(f, collection.items)
    # Rasmus: Collects the filtered items into a new SlotCollection
    SlotCollection(collect(filtered_items))
end

# Filtering functions
function get_slots_by_deck(collection::SlotCollection, deck_id::Int)
    filter(slot -> slot.deck_id == deck_id, collection.items)
end

function get_slots_by_cargo_type(collection::SlotCollection, cargo_type_id::Int)
    filter(slot -> slot.cargo_type_id == cargo_type_id, collection.items)
end

function get_refrigerated_slots(collection::SlotCollection)
    filter(slot -> slot.refrigerated, collection.items)
end

function get_deck_slots(collection::SlotCollection)
    filter(slot -> slot.on_deck, collection.items)
end

# Pretty printing for collection
function Base.show(io::IO, ::MIME"text/plain", collection::SlotCollection)
    println(io, "SlotCollection:")
    println(io, "├─ Total slots: $(collection.total_slots)")
    println(io, "└─ Items:")
    for (i, slot) in enumerate(collection.items)
        prefix = i == length(collection.items) ? "   └─ " : "   ├─ "
        println(io, prefix, slot)
    end
end

# Test data generation
function generate_random_slot(id::Int, deck_id::Int)
    cargo_types = 1:5  # Matching the 5 cargo types from cargo.jl
    
    Slot(
        id=id,
        loadmaster_id=id,  # Assuming loadmaster_id is same as id for test data
        deck_id=deck_id,
        cargo_type_id=rand(cargo_types),
        lcg=rand(0.0:0.1:100.0),
        tcg=rand(-10.0:0.1:10.0),
        vcg=rand(0.0:0.1:20.0),
        refrigerated=rand(Bool),
        length=rand([4.5, 12.0, 20.0, 40.0, 10.0]),  # Matching cargo dimensions
        width=rand([1.8, 2.5, 2.4, 2.4, 2.5]),      # Matching cargo dimensions
        on_deck=rand(Bool)
    )
end

function generate_test_slots(num_slots::Int, num_decks::Int)
    slots = Slot[]
    for deck_id in 1:num_decks
        slots_per_deck = num_slots ÷ num_decks
        for i in 1:slots_per_deck
            push!(slots, generate_random_slot(length(slots) + 1, deck_id))
        end
    end
    SlotCollection(slots)
end


function get_cargo_type_id(filename::String)
    # Extract cargo type from filename (without extension)
    cargo_name = splitext(basename(filename))[1]
    cargo_types = Dict(
        "Trailer" => 1,
        "Cars" => 2,
        "Machine" => 3,
        "Secu" => 4,
    )
    return get(cargo_types, cargo_name, 1)  # Default to 1 if not found
end

function get_deck_id(deck_name::String)
    deck_types = Dict(
        "UDECK" => 1,
        "MDECK" => 2,
        "TTOP" => 3
    )
    return get(deck_types, deck_name, 1)  # Default to 1 if not found
end

# Rasmus: Loads slots from one excel file.
# There is one excel file for each type of cargo.
function load_slots_from_excel(filepath::String, start_id::Int)
    xf = XLSX.readxlsx(filepath)
    sheet = xf[2]  # Use first sheet
    cargo_type_id = get_cargo_type_id(filepath)
    
    slots = Slot[]
    current_id = start_id
    
    for row in XLSX.eachrow(sheet)
        if row[1] == "G_RefNo" continue end
        
        slot = Slot(
            id = current_id,
            loadmaster_id = row[1],
            deck_id = get_deck_id(row[15]),
            cargo_type_id = cargo_type_id,  # Use cargo type from filename
            lcg = row[6],
            tcg = row[7],
            vcg = row[8],
            refrigerated = false,
            length = parse(Float64, replace(string(row[4]), "," => ".")),
            width = parse(Float64, replace(string(row[5]), "," => ".")),
            on_deck = row[15] == "UDECK"
        )
        push!(slots, slot)
        current_id += 1
    end
    
    return slots
end

# Rasmus: Loads all slots from a directory
function load_all_slots(directory::String="data/slots/finlandia")
    all_slots = Slot[]
    current_id = 1
    
    # Check if directory exists
    if !isdir(directory)
        error("Directory not found: $directory")
    end
    
    # Process all Excel files
    for file in readdir(directory)
        if endswith(lowercase(file), ".xls") || endswith(lowercase(file), ".xlsx")
            filepath = joinpath(directory, file)
            try
                slots = load_slots_from_excel(filepath, current_id)
                append!(all_slots, slots)
                current_id += length(slots)
            catch e
                @warn "Error processing $file: $e"
            end
        end
    end
    
    return SlotCollection(all_slots)
end