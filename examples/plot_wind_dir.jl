# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Status:
# Works fine for one to nine wind turbines, not clear how to plot more wind turbines.

using ControlPlots, FLORIDyn

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=true, save=true, rel_v_min=20.0, up_int = 4)

@assert Threads.nthreads() == 1 "This script is written for single threaded operation."

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic threading/parallel detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Arrays to store time series data
times = Float64[]
wind_directions = Vector{Float64}[]
turbines = 1:min(9, wf.nT)  # Use actual number of turbines, max 9

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
plot_data = [wind_dir_matrix[:, i] for i in 1:n_turbines]
turbine_labels = ["dir_t$i [°]" for i in turbines]

# Plot with dynamic number of turbines
p = plotx(times, plot_data...; fig="Wind Direction", xlabel="rel_time [s]", ylabels=turbine_labels, ysize = 10)
display(p)

nothing