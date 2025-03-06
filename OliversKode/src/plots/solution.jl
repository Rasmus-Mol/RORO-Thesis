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