# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Calculate axial induction factor, and calculate the demand

using FLORIDyn, ControlPlots, YAML

settings_file = get_default_project()[2]

# get the settings for the wind field, simulator, controller and turbine array
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# Calculate time step and relative end time based on simulation settings
time_step = sim.time_step  # seconds
t_end = sim.end_time - sim.start_time  # relative end time in seconds

function calc_axial_induction(turbine_group, time)
    # Example: simple constant axial induction factor
    return 1/3  # Constant value for all turbines
end

function calc_axial_induction(turbine, time)
    turbine_group = turbine_group(ta, turbine)
    return 1/3  # Constant value for all turbines
end

function calc_demand(time)
    # Example: time-varying demand with sinusoidal pattern
    initial_demand = 0.5
    final_demand = 1.0
    t1 = 240.0  # Time to start increasing demand
    t2 = 960.0  # Time to reach final demand
    if time < t1
        return initial_demand
    elseif time < t2
        return initial_demand + (final_demand - initial_demand) * (time - t1) / (t2 - t1)
    else
        return final_demand
    end
end

function plot_demand()
    # Create time vector from 0 to t_end with time_step intervals
    time_vector = 0:time_step:t_end
    
    # Calculate demand for each time point
    demand_values = [calc_demand(t) for t in time_vector]
    
    # Create the plot
    ControlPlots.plot(time_vector, demand_values, 
                     xlabel="Time [s]", 
                     ylabel="Demand [-]", 
                     title="Demand vs Time")
end