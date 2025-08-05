# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Status:
# Works fine for one to nine wind turbines, not clear how to plot more wind turbines.

using ControlPlots, FLORIDyn

# Dialog to set MULTI variable
println("\033[1mPlot wind direction for multiple turbines?\033[0m")
print("Enter 'y' for multiple turbines, 'n' for single turbine [y/N]: ")
response = readline()
MULTI = lowercase(strip(response)) in ["y", "yes", "true", "1"]

if MULTI
    println("Multi-turbine mode selected, up to 81 turbines.")
else
    println("Single turbine per subplot selected, up to 9 turbines.")
end

settings_file = "data/2021_9T_Data.yaml"

@assert Threads.nthreads() == 1 "This script is written for single threaded operation.
                                  Quit Julia and start it with 'jl2'."

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic threading/parallel detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Arrays to store time series data
times = Float64[]
wind_directions = Vector{Float64}[]
turbines = 1:wf.nT

for time in sim.start_time:sim.time_step:sim.end_time
    local wind_direction
    rel_time = time - sim.start_time
    
    # Calculate wind direction at current time
    # Try to get time-varying wind direction using the direction model
    wind_dir_vec = getWindDirT(set.dir_mode, wind.dir, turbines, time)
    wind_direction = wind_dir_vec  # Vector for all requested turbines
    
    # Store the data
    push!(times, rel_time)
    push!(wind_directions, wind_direction)
end

# Convert vector of vectors to matrix for easier plotting
wind_dir_matrix = hcat(wind_directions...)  # Transpose to get time × turbines
wind_dir_matrix = wind_dir_matrix'  # Now it's time × turbines

# Create dynamic plot arguments based on number of turbines
n_turbines = length(turbines)

if MULTI
    rows, lines = get_layout(wf.nT)
    
    # Group turbines into subplots based on layout
    plot_data = []
    turbine_labels = []
    subplot_labels = []
    
    local turbine_idx = 1
    for row in 1:rows
        if turbine_idx > n_turbines
            break
        end
        
        # Collect lines for this subplot
        local lines_in_subplot = Vector{Vector{Float64}}()
        local labels_in_subplot = Vector{String}()
        
        for line in 1:lines
            if turbine_idx <= n_turbines
                push!(lines_in_subplot, wind_dir_matrix[:, turbine_idx])
                push!(labels_in_subplot, "T$(turbines[turbine_idx])")
                turbine_idx += 1
            end
        end
        
        # Add subplot data
        if length(lines_in_subplot) == 1
            push!(plot_data, lines_in_subplot[1])
        else
            push!(plot_data, lines_in_subplot)
        end
        
        push!(turbine_labels, "Wind Direction [°]")
        push!(subplot_labels, labels_in_subplot)
    end
    
    # Plot with multiple lines per subplot
    p = plotx(times, plot_data...; ylabels=turbine_labels, labels=subplot_labels,
              fig="Wind Direction", xlabel="rel_time [s]", ysize = 10, bottom=0.02)
else
    # Single turbine mode - one turbine per subplot
    plot_data = [wind_dir_matrix[:, i] for i in 1:n_turbines]
    turbine_labels = ["T$i wind_dir [°]" for i in turbines]
    
    p = plotx(times, plot_data...; fig="Wind Direction", xlabel="rel_time [s]", 
              ylabels=turbine_labels, ysize = 9, bottom=0.02)
end

display(p)

nothing