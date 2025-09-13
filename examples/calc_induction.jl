# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Calculate axial induction factor, and calculate the demand

# TODO
# - fill the field con.induction_data with a suitable Matrix
# - use this induction data in the the function get_axial_induction

using FLORIDyn, ControlPlots, YAML

const cp_max = 16/27  # Betz limit

settings_file = get_default_project()[2]

# get the settings for the wind field, simulator, controller and turbine array
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# Calculate time step and relative end time based on simulation settings
time_step = sim.time_step  # seconds
t_end = sim.end_time - sim.start_time  # relative end time in seconds

function calc_induction_per_group(turbine_group, time)

    # simple example: assume no wakes
    demand = calc_demand(time)
    induction = calc_induction(demand * cp_max)
    return induction
end

function calc_axial_induction(ta, con, turbine, time)
    turbine_group = turbine_group(ta, turbine)
    return calc_induction_per_group(turbine_group, time)
end

function calc_cp(induction)
    return 4 * induction * (1 - induction)^2
end
function calc_induction(cp)
    # Solve the equation: cp = 4*a*(1-a)^2 for induction factor 'a'
    
    # Check if cp is within valid range
    if cp < 0 || cp > cp_max
        error("cp value $cp is outside valid range [0, $(cp_max)]")
    end
    
    # Special cases
    if cp ≈ 0
        return 0.0
    elseif cp ≈ cp_max
        return 1/3
    end
    
    # Use a simple bisection method which is more robust
    f(a) = 4*a*(1-a)^2 - cp
    
    # Find the interval containing the root
    # We know the function has roots in [0, 1/3] and [1/3, 1] for most cp values
    a_low, a_high = 0.0, 1.0
    
    # Check if we need to search in the lower or upper interval
    if f(1/3) >= 0
        # Root is in [0, 1/3]
        a_high = 1/3
    else
        # Root is in [1/3, 1]
        a_low = 1/3
    end
    
    # Bisection method
    tolerance = 1e-12
    max_iterations = 100
    
    for i in 1:max_iterations
        a_mid = (a_low + a_high) / 2
        f_mid = f(a_mid)
        
        if abs(f_mid) < tolerance || (a_high - a_low) < tolerance
            return a_mid
        end
        
        if f(a_low) * f_mid < 0
            a_high = a_mid
        else
            a_low = a_mid
        end
        
        if i == max_iterations
            error("Bisection method did not converge for cp=$cp")
        end
    end
    
    return (a_low + a_high) / 2
end

function calc_demand(time)
    # Example: linearly increasing demand from 0.5 to 1.0 over the simulation time
    initial_demand = 0.5
    final_demand = 0.8
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
    induction_values = [calc_induction_per_group(1, t) for t in time_vector]  # Example for group 1
    
    # Create the plot
    ControlPlots.plot(time_vector, [demand_values, induction_values], 
                     xlabel="Time [s]", 
                     ylabel="Demand and Induction [-]", 
                     title="Demand and Induction vs Time")
end