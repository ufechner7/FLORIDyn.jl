# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Plot turbine positions from windFarmProperties.yaml

using ControlPlots, YAML

yaml_file = joinpath(@__DIR__, "..", "data", "windFarmProperties.yaml")

function plot_wind_farm(yaml_file; n_groups=9)
    # Load YAML data
    data = YAML.load_file(yaml_file)
    turbines = data["turbines"]

    # Extract turbine names and positions
    turbine_names = sort(collect(keys(turbines)))  # Sorted: T1, T10, T100, ...
    # Sort naturally by turbine number
    turbine_names = sort(turbine_names; by = t -> parse(Int, t[2:end]))

    n_turbines = length(turbine_names)
    x_coords = zeros(n_turbines)
    y_coords = zeros(n_turbines)

    for (i, name) in enumerate(turbine_names)
        pos = turbines[name]["baseLocation"]
        x_coords[i] = pos[1]
        y_coords[i] = pos[2]
    end

    # If n_groups > 0, assign groups based on X coordinate
    if n_groups > 0
        sorted_x = sort(x_coords)
        group_edges = [sorted_x[1]]
        step = n_turbines ÷ n_groups
        for g in 1:(n_groups - 1)
            push!(group_edges, sorted_x[g * step])
        end
        push!(group_edges, sorted_x[end] + 1)

        group_ids = zeros(Int, n_turbines)
        for i in 1:n_turbines
            for g in 1:n_groups
                if x_coords[i] >= group_edges[g] && x_coords[i] < group_edges[g + 1]
                    group_ids[i] = g
                    break
                end
            end
        end
    else
        group_ids = ones(Int, n_turbines)
        n_groups = 1
    end

    # Colors and markers for groups
    colors = ["red", "blue", "green", "orange", "purple", "brown", "pink", "gray", "olive"]
    markers = ["^", "o", "s", "D", "v", "<", ">", "p", "*"]

    # Create plot
    plt.figure(figsize=(15, 6))

    for g in 1:n_groups
        mask = group_ids .== g
        gx = x_coords[mask]
        gy = y_coords[mask]
        color = colors[mod1(g, length(colors))]
        marker = markers[mod1(g, length(markers))]
        plt.scatter(gx, gy, s=120, c=color, marker=marker,
                   label="TG$g ($(count(mask)) turbines)",
                   alpha=0.8, edgecolors="black", linewidth=0.5)
    end

    # Add turbine labels
    for i in 1:n_turbines
        plt.annotate(turbine_names[i], (x_coords[i], y_coords[i]), xytext=(5, 5),
                    textcoords="offset points", fontsize=8,
                    bbox=Dict("boxstyle"=>"round,pad=0.2", "facecolor"=>"white", "alpha"=>0.7))
    end

    plt.xlabel("X Coordinate [m]")
    plt.ylabel("Y Coordinate [m]")
    plt.title("Wind Farm Layout - $n_groups Turbine Groups by X Coordinate")
    plt.grid(true, alpha=0.3)
    plt.legend(loc="upper center", bbox_to_anchor=(0.69, 0.4))

    # Margins
    x_range = maximum(x_coords) - minimum(x_coords)
    y_range = maximum(y_coords) - minimum(y_coords)
    x_margin = x_range * 0.05
    y_margin = y_range * 0.05
    plt.xlim(minimum(x_coords) - x_margin, maximum(x_coords) + x_margin)
    plt.ylim(minimum(y_coords) - y_margin, maximum(y_coords) + y_margin)

    plt.tight_layout(pad=0.5)
    plt.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.92)

    # Print statistics
    println("Wind Farm Statistics:")
    println("  Number of turbines: $n_turbines")
    println("  X range: $(round(minimum(x_coords), digits=1)) - $(round(maximum(x_coords), digits=1)) m")
    println("  Y range: $(round(minimum(y_coords), digits=1)) - $(round(maximum(y_coords), digits=1)) m")
    println("  Farm width: $(round(x_range, digits=1)) m")
    println("  Farm height: $(round(y_range, digits=1)) m")

    return plt.gcf()
end

# Run the plotting function
println("Loading turbine data from: $yaml_file")
try
    fig = plot_wind_farm(yaml_file; n_groups=9)
    println("\nPlot created successfully with group colors!")
catch e
    println("Error loading or plotting turbine data: $e")
    rethrow()
end