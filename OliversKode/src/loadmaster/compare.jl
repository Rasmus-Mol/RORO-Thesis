# Define struct for LoadMaster results
@kwdef struct LoadMasterResults
    displacement::Float64
    dwt::Float64
    draft_fp::Float64
    draft_mid::Float64
    draft_ap::Float64
    heel::Float64
    trim::Float64
    gm_actual::Float64
    gm_min::Float64
end

# Parse LoadMaster items into struct
function parse_loadmaster_results(items::Vector{Item})
    # Helper to find value by description
    function find_value(items, desc)
        item = first(filter(i -> i.description == desc, items))
        parse(Float64, replace(item.value, ['p', 'f'] .=> ""))
    end

    LoadMasterResults(
        displacement = find_value(items, "Displ"),
        dwt = find_value(items, "DWT"),
        draft_fp = find_value(items, "dfp"),
        draft_mid = find_value(items, "dm"),
        draft_ap = find_value(items, "dap"),
        heel = find_value(items, "Heel"),
        trim = find_value(items, "Trim"),
        gm_actual = find_value(items, "G'Mact0"),
        gm_min = find_value(items, "G'Mmin")
    )
end

# Compare solution with LoadMaster results
# Rasmus: I think this should be Solution1 instead of Solution1
#function compare_solution(sol::Solution1, lm_results::LoadMasterResults)
function compare_solution(sol::Solution, lm_results::LoadMasterResults)
    # Calculate differences
    weight_diff = abs(sol.total_weight - lm_results.displacement)
    lcg_impact = abs(sol.lcg - lm_results.trim)
    tcg_impact = abs(sol.tcg - lm_results.heel)

    # Return comparison results
    (
        weight_difference = weight_diff,
        weight_match = weight_diff < 1.0,  # Within 1 ton
        trim_difference = lcg_impact,
        trim_match = lcg_impact < 0.1,     # Within 0.1m
        heel_difference = tcg_impact,
        heel_match = tcg_impact < 0.1      # Within 0.1 degrees
    )
end

# Usage example:
lm_results = parse_loadmaster_results(result)
comparison = compare_solution(solution, lm_results)