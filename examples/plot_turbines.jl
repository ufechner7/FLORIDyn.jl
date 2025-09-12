# plot the turbines of the given project

using FLORIDyn, ControlPlots, YAML

settings_file = get_default_project()[2]

function plot_turbines(ta::TurbineArray)
    # Extract x and y coordinates from the position matrix
    x_coords = ta.pos[:, 1]  # First column contains x coordinates
    y_coords = ta.pos[:, 2]  # Second column contains y coordinates
    
    # Create scatter plot of turbine locations
    plt.figure(figsize=(10, 8))
    plt.scatter(x_coords, y_coords, s=100, c="red", marker="^", label="Wind Turbines", alpha=0.7)
    
    # Add turbine IDs as labels
    for i in 1:length(x_coords)
        plt.annotate("T$i", (x_coords[i], y_coords[i]), xytext=(5, 5), 
                    textcoords="offset points", fontsize=8)
    end
    
    plt.xlabel("X Coordinate [m]")
    plt.ylabel("Y Coordinate [m]")
    plt.title("Wind Farm Layout - Turbine Locations")
    plt.grid(true, alpha=0.3)
    plt.legend()
    plt.axis("equal")
    plt.tight_layout()
    
    return plt.gcf()
end

# Load the settings and extract turbine data
println("Loading turbine data from: $settings_file")
try
    # Use setup function to load all configuration data including turbine array
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    
    println("Found $(size(ta.pos, 1)) turbines in the configuration")
    println("Turbine types: $(unique(ta.type))")
    
    # Plot the turbines
    fig = plot_turbines(ta)
    
    # Display some statistics
    x_coords = ta.pos[:, 1]
    y_coords = ta.pos[:, 2]
    println("\nWind Farm Statistics:")
    println("  X range: $(round(minimum(x_coords), digits=1)) - $(round(maximum(x_coords), digits=1)) m")
    println("  Y range: $(round(minimum(y_coords), digits=1)) - $(round(maximum(y_coords), digits=1)) m")
    println("  Farm width: $(round(maximum(x_coords) - minimum(x_coords), digits=1)) m")
    println("  Farm height: $(round(maximum(y_coords) - minimum(y_coords), digits=1)) m")
    
    println("\nPlot created successfully!")
    
catch e
    println("Error loading or plotting turbine data: $e")
    println("Please check that the settings file exists and contains valid turbine data.")
end