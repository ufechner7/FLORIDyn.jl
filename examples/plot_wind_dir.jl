# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using DistributedNext, Timers, ControlPlots, FLORIDyn

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=true, save=true, rel_v_min=20.0, up_int = 4)

# Automatic parallel/threading setup
include("../src/visualisation/smart_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic threading/parallel detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Arrays to store time series data
times = Float64[]
wind_speeds = Float64[]

for time in sim.start_time:sim.time_step:sim.end_time
    local wind_speed
    rel_time = time - sim.start_time
    
    # Calculate wind speed at current time
    # Use the wind velocity function with the current time
    try
        # Try to get time-varying wind speed using the velocity model
        wind_speed_vec = getWindSpeedT(set.vel_mode, wind.vel, [1], time)
        wind_speed = wind_speed_vec[1]  # Extract first value since we requested turbine 1
    catch e
        @warn "Failed to get time-varying wind speed: $e"
        # Fallback to constant wind speed
        if hasfield(typeof(wind.vel), :Data)
            wind_speed = wind.vel.Data
        else
            wind_speed = wind.vel  # For constant velocity
        end
    end
    
    # Store the data
    push!(times, rel_time)
    push!(wind_speeds, wind_speed)
    
    println("Time: $(rel_time)s, Wind Speed: $(round(wind_speed, digits=2)) m/s")
end
# plot(times, wind_speeds)

nothing