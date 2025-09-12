# plot the turbines of the given project

using FLORIDyn, ControlPlots, YAML

settings_file = get_default_project()[2]

GROUPS = 4

function plot_turbines(ta::TurbineArray, turbine_groups)
    # Extract x and y coordinates from the position matrix
    x_coords = ta.pos[:, 1]  # First column contains x coordinates
    y_coords = ta.pos[:, 2]  # Second column contains y coordinates
    
    # Define colors for each group
    colors = ["red", "blue", "green", "orange", "purple", "brown", "pink", "gray"]
    markers = ["^", "o", "s", "D", "v", "<", ">", "p"]
    
    # Create scatter plot of turbine locations
    plt.figure(figsize=(15, 6))  # Wider, shorter figure to match data proportions
    
    # Plot each group with different colors
    plotted_groups = Set{String}()
    
    for group in turbine_groups
        if group["name"] == "all"
            continue  # Skip the "all" group
        end
        
        group_turbines = group["turbines"]
        group_name = group["name"]
        group_id = group["id"]
        
        # Get coordinates for this group
        group_x = [x_coords[t] for t in group_turbines]
        group_y = [y_coords[t] for t in group_turbines]
        
        # Select color and marker for this group
        color = colors[min(group_id, length(colors))]
        marker = markers[min(group_id, length(markers))]
        
        # Plot this group
        plt.scatter(group_x, group_y, s=120, c=color, marker=marker, 
                   label="$group_name ($(length(group_turbines)) turbines)", 
                   alpha=0.8, edgecolors="black", linewidth=0.5)
        
        push!(plotted_groups, group_name)
    end
    
    # Add turbine IDs as labels
    for i in 1:length(x_coords)
        plt.annotate("T$i", (x_coords[i], y_coords[i]), xytext=(5, 5), 
                    textcoords="offset points", fontsize=8, 
                    bbox=Dict("boxstyle"=>"round,pad=0.2", "facecolor"=>"white", "alpha"=>0.7))
    end
    
    plt.xlabel("X Coordinate [m]")
    plt.ylabel("Y Coordinate [m]")
    plt.title("Wind Farm Layout - Turbine Groups by X Coordinate")
    plt.grid(true, alpha=0.3)
    plt.legend(loc="lower right")
    
    # Set axis limits to reduce unused space
    y_min = minimum(y_coords)
    y_max = maximum(y_coords)
    y_range = y_max - y_min
    y_margin = y_range * 0.05  # 5% margin on each side
    
    x_min = minimum(x_coords)
    x_max = maximum(x_coords)
    x_range = x_max - x_min
    x_margin = x_range * 0.05  # 5% margin on each side

    plt.xlim(x_min - x_margin, x_max + x_margin)
    plt.ylim(y_min - y_margin, y_max + y_margin)
    
    # Don't force equal aspect ratio to allow better use of figure space
    plt.tight_layout(pad=0.5)
    
    # Further reduce margins around the plot
    plt.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.92)
    
    return plt.gcf()
end

# Load the settings and extract turbine data
println("Loading turbine data from: $settings_file")
try
    # Load YAML data to get turbine groups
    yaml_data = YAML.load_file(settings_file)
    turbine_groups = yaml_data["turbine_groups"]
    
    # Use setup function to load all configuration data including turbine array
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    
    println("Found $(size(ta.pos, 1)) turbines in the configuration")
    println("Turbine types: $(unique(ta.type))")
    
    # Display group information
    println("\nTurbine Groups:")
    for group in turbine_groups
        if group["name"] != "all"
            println("  $(group["name"]): $(length(group["turbines"])) turbines")
        end
    end
    
    # Plot the turbines with group colors
    fig = plot_turbines(ta, turbine_groups)
    
    # Display some statistics
    x_coords = ta.pos[:, 1]
    y_coords = ta.pos[:, 2]
    println("\nWind Farm Statistics:")
    println("  X range: $(round(minimum(x_coords), digits=1)) - $(round(maximum(x_coords), digits=1)) m")
    println("  Y range: $(round(minimum(y_coords), digits=1)) - $(round(maximum(y_coords), digits=1)) m")
    println("  Farm width: $(round(maximum(x_coords) - minimum(x_coords), digits=1)) m")
    println("  Farm height: $(round(maximum(y_coords) - minimum(y_coords), digits=1)) m")
    
    println("\nPlot created successfully with group colors!")
    
catch e
    println("Error loading or plotting turbine data: $e")
    println("Please check that the settings file exists and contains valid turbine data.")
end