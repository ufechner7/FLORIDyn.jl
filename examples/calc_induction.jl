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

function calc_cp(induction)
    return 4 * induction * (1 - induction)^2
end
function calc_induction(cp)
    # Solve cubic equation: 4a^3 - 8a^2 + 4a - cp = 0
    # Using Cardano's method for cubic equations
    a = 0.0
    b = -2.0
    c = 1.0
    d = -cp / 4.0

    Δ0 = b^2 - 3*a*c
    Δ1 = 2*b^3 - 9*a*b*c + 27*a^2*d

    C = ((Δ1 + sqrt(Δ1^2 - 4*Δ0^3)) / 2)^(1/3)

    if C == 0
        C = ((Δ1 - sqrt(Δ1^2 - 4*Δ0^3)) / 2)^(1/3)
    end

    ξ = exp(2π*im/3)  # Complex cube root of unity

    roots = [
        -(1/(3*a)) * (b + C + Δ0 / C),
        -(1/(3*a)) * (b + ξ*C + Δ0 / (ξ*C)),
        -(1/(3*a)) * (b + ξ^2*C + Δ0 / (ξ^2*C))
    ]

    # Filter real roots in the range [0, 1]
    real_roots = filter(r -> isreal(r) && r >= 0 && r <= 1, roots)
    
    if length(real_roots) == 0
        error("No valid real root found for cp=$cp")
    end

    return real(real_roots[1])  # Return the first valid real root
end

function calc_demand(time)
    # Example: linearly increasing demand from 0.5 to 1.0 over the simulation time
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