function plot_solution(solution::Solution)
    n_decks = maximum(cargo.deck for cargo in solution.cargo)
    plots = []

    min_weight = minimum(c.weight for c in solution.cargo)
    max_weight = maximum(c.weight for c in solution.cargo)

    colorgrad = cgrad(:roma, rev = true)

    # Calculate bounds based on vessel dimensions
    deck_bounds = map(1:n_decks) do deck
        deck_cargo = filter(c -> c.deck == deck, solution.cargo)
        if isempty(deck_cargo)
            return (0.0, 1.0, 0.0, 1.0)  # Default bounds for empty decks
        end
        min_x = minimum(c.lcg - c.length/2 for c in deck_cargo)
        max_x = maximum(c.lcg + c.length/2 for c in deck_cargo)
        min_y = minimum(c.tcg - c.width/2 for c in deck_cargo)
        max_y = maximum(c.tcg + c.width/2 for c in deck_cargo)
        return (min_x, max_x, min_y, max_y)
    end

    # Global bounds for consistent scaling
    global_min_x = minimum(b[1] for b in deck_bounds)
    global_max_x = maximum(b[2] for b in deck_bounds)
    global_min_y = minimum(b[3] for b in deck_bounds)
    global_max_y = maximum(b[4] for b in deck_bounds)

    margin = 5.0

    for deck in 1:n_decks
        # Create empty plot with correct bounds
        p = plot(
            xlim = (global_min_x - margin, global_max_x + margin),
            ylim = (global_min_y - margin, global_max_y + margin),
            xlabel = "Deck $deck - Length (m)",
            ylabel = "Width (m)",
            grid = true,
            aspect_ratio = :equal,
            right_margin = 10Plots.mm,
        )
        if deck == 1
            title!(p, "Cargo Placement - Deterministic Model")
        end
        # Add cargo boxes for this deck
        deck_cargo = filter(c -> c.deck == deck, solution.cargo)

        for cargo in deck_cargo
            x = cargo.lcg - cargo.length/2
            y = cargo.tcg - cargo.width/2
            w = cargo.width
            l = cargo.length

            rect = Shape([
                (x, y),
                (x + l, y),
                (x + l, y + w),
                (x, y + w),
            ])

            plot!(p, rect,
                fillcolor = get(colorgrad, cargo.weight, (min_weight, max_weight)),
                alpha = 0.7,
                label = false,
                linewidth = 1,
            )

            # Add cargo ID label
            annotate!(p,
                x + l/2,
                y + w/2,
                text("$(cargo.id)", 10, :center),
            )
        end

        push!(plots, p)
    end

    # plots = reverse(plots)  # Reverse order for correct deck numbering

    # Add colorbar
    colorbar = scatter([0, 0], [0, 1], zcolor = [0, 3],
        clims = (min_weight, max_weight),
        xlims = (1, 1.1),
        label = "",
        c = colorgrad,
        colorbar_title = "Cargo Weight",
        framestyle = :none
    )
    push!(plots, colorbar)

    l = @layout [grid(n_decks, 1) a{0.001w}]

    final_plot = plot(plots...,
        size = (1500, 300 * n_decks),
        layout = l,
        link = :all,
    )
    
    return final_plot
end

function plot_solution_random_plan(solution::Solution)
    n_decks = maximum(cargo.deck for cargo in solution.cargo)
    plots = []

    min_weight = minimum(c.weight for c in solution.cargo)
    max_weight = maximum(c.weight for c in solution.cargo)

    colorgrad = cgrad(:roma, rev = true)

    # Calculate bounds based on vessel dimensions
    deck_bounds = map(1:n_decks) do deck
        deck_cargo = filter(c -> c.deck == deck, solution.cargo)
        if isempty(deck_cargo)
            return (0.0, 1.0, 0.0, 1.0)  # Default bounds for empty decks
        end
        min_x = minimum(c.lcg - c.length/2 for c in deck_cargo)
        max_x = maximum(c.lcg + c.length/2 for c in deck_cargo)
        min_y = minimum(c.tcg - c.width/2 for c in deck_cargo)
        max_y = maximum(c.tcg + c.width/2 for c in deck_cargo)
        return (min_x, max_x, min_y, max_y)
    end

    # Global bounds for consistent scaling
    global_min_x = minimum(b[1] for b in deck_bounds)
    global_max_x = maximum(b[2] for b in deck_bounds)
    global_min_y = minimum(b[3] for b in deck_bounds)
    global_max_y = maximum(b[4] for b in deck_bounds)

    margin = 5.0

    for deck in 1:n_decks
        # Create empty plot with correct bounds
        p = plot(
            xlim = (global_min_x - margin, global_max_x + margin),
            ylim = (global_min_y - margin, global_max_y + margin),
            xlabel = "Deck $deck - Length (m)",
            ylabel = "Width (m)",
            grid = true,
            aspect_ratio = :equal,
            right_margin = 10Plots.mm,
        )
        if deck == 1
            title!(p, "Cargo Placement - Random stowage plan")
        end
        # Add cargo boxes for this deck
        deck_cargo = filter(c -> c.deck == deck, solution.cargo)

        for cargo in deck_cargo
            x = cargo.lcg - cargo.length/2
            y = cargo.tcg - cargo.width/2
            w = cargo.width
            l = cargo.length

            rect = Shape([
                (x, y),
                (x + l, y),
                (x + l, y + w),
                (x, y + w),
            ])

            plot!(p, rect,
                fillcolor = get(colorgrad, cargo.weight, (min_weight, max_weight)),
                alpha = 0.7,
                label = false,
                linewidth = 1,
            )

            # Add cargo ID label
            annotate!(p,
                x + l/2,
                y + w/2,
                text("$(cargo.id)", 10, :center),
            )
        end

        push!(plots, p)
    end

    # plots = reverse(plots)  # Reverse order for correct deck numbering

    # Add colorbar
    colorbar = scatter([0, 0], [0, 1], zcolor = [0, 3],
        clims = (min_weight, max_weight),
        xlims = (1, 1.1),
        label = "",
        c = colorgrad,
        colorbar_title = "Cargo Weight",
        framestyle = :none
    )
    push!(plots, colorbar)

    l = @layout [grid(n_decks, 1) a{0.001w}]

    final_plot = plot(plots...,
        size = (1500, 300 * n_decks),
        layout = l,
        link = :all,
    )
    
    return final_plot
end

# Plot a realization of a scenario
function plot_solution_stochastic(solution::SolutionStochastic,sc)
    n_decks = maximum(cargo.deck for cargo in  solution.cargo)
    plots = []
    # Get the weights of the cargos in the scenario
    cargo_weights = []
    for i in 1:solution.n_cargo_loaded
        cargo_weights = push!(cargo_weights, solution.cargo[i].weight[sc])
    end
    min_weight = minimum(cargo_weights)
    max_weight = maximum(cargo_weights)

    colorgrad = cgrad(:roma, rev = true)

    # Calculate bounds based on vessel dimensions
    deck_bounds = map(1:n_decks) do deck
        deck_cargo = filter(c -> c.deck == deck, solution.cargo)
        if isempty(deck_cargo)
            return (0.0, 1.0, 0.0, 1.0)  # Default bounds for empty decks
        end
        min_x = minimum(c.lcg - c.length/2 for c in deck_cargo)
        max_x = maximum(c.lcg + c.length/2 for c in deck_cargo)
        min_y = minimum(c.tcg - c.width/2 for c in deck_cargo)
        max_y = maximum(c.tcg + c.width/2 for c in deck_cargo)
        return (min_x, max_x, min_y, max_y)
    end

    # Global bounds for consistent scaling
    global_min_x = minimum(b[1] for b in deck_bounds)
    global_max_x = maximum(b[2] for b in deck_bounds)
    global_min_y = minimum(b[3] for b in deck_bounds)
    global_max_y = maximum(b[4] for b in deck_bounds)

    margin = 5.0

    for deck in 1:n_decks
        # Create empty plot with correct bounds
        p = plot(
            xlim = (global_min_x - margin, global_max_x + margin),
            ylim = (global_min_y - margin, global_max_y + margin),
            xlabel = "Deck $deck - Length (m)",
            ylabel = "Width (m)",
            grid = true,
            aspect_ratio = :equal,
            right_margin = 10Plots.mm,
        )
        if deck == 1
            title!(p, "Cargo Placement - Stochastic Model, scenario $sc")
        end
        # Add cargo boxes for this deck
        deck_cargo = filter(c -> c.deck == deck, solution.cargo)

        for cargo in deck_cargo
            x = cargo.lcg - cargo.length/2
            y = cargo.tcg - cargo.width/2
            w = cargo.width
            l = cargo.length

            rect = Shape([
                (x, y),
                (x + l, y),
                (x + l, y + w),
                (x, y + w),
            ])

            plot!(p, rect,
                fillcolor = get(colorgrad, cargo.weight[sc], (min_weight, max_weight)),
                alpha = 0.7,
                label = false,
                linewidth = 1,
            )

            # Add cargo ID label
            annotate!(p,
                x + l/2,
                y + w/2,
                text("$(cargo.id)", 10, :center),
            )
        end

        push!(plots, p)
    end

    # plots = reverse(plots)  # Reverse order for correct deck numbering

    # Add colorbar
    colorbar = scatter([0, 0], [0, 1], zcolor = [0, 3],
        clims = (min_weight, max_weight),
        xlims = (1, 1.1),
        label = "",
        c = colorgrad,
        colorbar_title = "Cargo Weight",
        framestyle = :none
    )
    push!(plots, colorbar)

    l = @layout [grid(n_decks, 1) a{0.001w}]

    final_plot = plot(plots...,
        size = (1500, 300 * n_decks),
        layout = l,
        link = :all,
    )
    return final_plot
end

# Plot ballast weight difference between stochastic and deterministic model
function plot_ballast_weight_diff(diff_ballast_weight,sc,n)
    p = plot(
        diff_ballast_weight,
        xlabel = "Scenario",
        ylabel = "Difference in ballast weight (t)",
        title = "Difference in ballast weight\nbetween stochastic and deterministic model,\n scenarios = $sc, unknown weight cargo = $n",
        legend = false
    )
    return p
end
# Plot ballast and cargo weight difference between stochastic and deterministic model
function plot_ballast_cargo_weight_diff(diff_ballast_weight,diff_cargo_weight,sc,n)
    p = plot(
        diff_ballast_weight,
        label = "Ballast",
        xlabel = "Scenario",
        ylabel = "Difference in weight (t)",
        title = "Difference in ballast and cargo weight\nbetween stochastic and deterministic model,\n scenarios = $sc, unknown weight cargo = $n",
        color =:blue
    )
    twinx()
    plot!(
        diff_cargo_weight,
        label = "Cargo",
        #ylabel = "Difference in cargo weight (t)",
        color =:red
    )
    return p
end

# Plots a random stowage plan.
function plot_solution_random_plan(cs, cargoC::CargoCollection,slots::SlotCollection)
    placements = CargoPlacement[]
    n_cargo_loaded = 0
    n_cargo = length(cs[:,1])
    n_slots = length(cs[1,:])
    for c in 1:n_cargo, s in 1:n_slots
        if cs[c,s] > 0.5  # Binary variable threshold
            n_cargo_loaded += 1
            slot = slots[s]
            # Wrong here
            cargo = filter(x -> x.id== c, cargoC)[1]
            #cargo = cargoC[c]
            push!(placements, CargoPlacement(
                id = cargo.id,
                cargo_type_id = cargo.cargo_type_id,
                deck = slot.deck_id,
                lcg = slot.lcg,
                tcg = slot.tcg,
                vcg = slot.vcg,
                weight = cargo.weight,
                length = get_length(cargo),
                width = get_width(cargo),
                height = get_height(cargo),
                haz_class = cargo.hazardous
            ))
        end
    end
    n_decks = maximum(cargo.deck for cargo in placements)
    plots = []

    min_weight = minimum(c.weight for c in cargoC)
    max_weight = maximum(c.weight for c in cargoC)

    colorgrad = cgrad(:roma, rev = true)

    # Calculate bounds based on vessel dimensions
    deck_bounds = map(1:n_decks) do deck
        deck_cargo = filter(c -> c.deck == deck, placements)
        if isempty(deck_cargo)
            return (0.0, 1.0, 0.0, 1.0)  # Default bounds for empty decks
        end
        min_x = minimum(c.lcg - c.length/2 for c in deck_cargo)
        max_x = maximum(c.lcg + c.length/2 for c in deck_cargo)
        min_y = minimum(c.tcg - c.width/2 for c in deck_cargo)
        max_y = maximum(c.tcg + c.width/2 for c in deck_cargo)
        return (min_x, max_x, min_y, max_y)
    end

    # Global bounds for consistent scaling
    global_min_x = minimum(b[1] for b in deck_bounds)
    global_max_x = maximum(b[2] for b in deck_bounds)
    global_min_y = minimum(b[3] for b in deck_bounds)
    global_max_y = maximum(b[4] for b in deck_bounds)

    margin = 5.0

    for deck in 1:n_decks
        # Create empty plot with correct bounds
        p = plot(
            xlim = (global_min_x - margin, global_max_x + margin),
            ylim = (global_min_y - margin, global_max_y + margin),
            xlabel = "Deck $deck - Length (m)",
            ylabel = "Width (m)",
            grid = true,
            aspect_ratio = :equal,
            right_margin = 10Plots.mm,
        )
        if deck == 1
            title!(p, "Cargo Placement - Random stowage plan")
        end
        # Add cargo boxes for this deck
        deck_cargo = filter(c -> c.deck == deck, placements)

        for cargo in deck_cargo
            x = cargo.lcg - cargo.length/2
            y = cargo.tcg - cargo.width/2
            w = cargo.width
            l = cargo.length

            rect = Shape([
                (x, y),
                (x + l, y),
                (x + l, y + w),
                (x, y + w),
            ])

            plot!(p, rect,
                fillcolor = get(colorgrad, cargo.weight, (min_weight, max_weight)),
                alpha = 0.7,
                label = false,
                linewidth = 1,
            )

            # Add cargo ID label
            annotate!(p,
                x + l/2,
                y + w/2,
                text("$(cargo.id)", 10, :center),
            )
        end

        push!(plots, p)
    end

    # plots = reverse(plots)  # Reverse order for correct deck numbering

    # Add colorbar
    colorbar = scatter([0, 0], [0, 1], zcolor = [0, 3],
        clims = (min_weight, max_weight),
        xlims = (1, 1.1),
        label = "",
        c = colorgrad,
        colorbar_title = "Cargo Weight",
        framestyle = :none
    )
    push!(plots, colorbar)

    l = @layout [grid(n_decks, 1) a{0.001w}]

    final_plot = plot(plots...,
        size = (1500, 300 * n_decks),
        layout = l,
        link = :all,
    )
    
    return final_plot
end
