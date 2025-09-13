# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Calculate axial induction factor, and calculate the demand

# TODO
# - integrate the function calc_axial_induction in FLORIDyn

using FLORIDyn, ControlPlots, YAML

const cp_max = 16/27  # Betz limit

settings_file = get_default_project()[2]

# get the settings for the wind field, simulator, controller and turbine array
wind, sim, con, floris, floridyn, ta = setup(settings_file)
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

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
    # Check if pre-calculated induction data is available
    if hasfield(typeof(con), :induction_data) && !isnothing(con.induction_data)
        # Use pre-calculated data from con.induction_data
        time_vector = 0:time_step:t_end
        
        # Find the time index that corresponds to the requested time
        time_idx = findfirst(t -> t >= time, time_vector)
        if isnothing(time_idx)
            time_idx = length(time_vector)  # Use last time step if time is beyond range
        end
        
        # Return the pre-calculated induction value for this turbine and time
        if turbine <= size(con.induction_data, 2) && time_idx <= size(con.induction_data, 1)
            return con.induction_data[time_idx, turbine]
        end
    end
    
    # Fallback to dynamic calculation if no pre-calculated data is available
    group_id = FLORIDyn.turbine_group(ta, turbine)
    
    # Apply corrections based on turbine group and time with linear interpolation
    # Rules:
    # - apply a large reduction (0.2) in the power for group 1 at t=0
    # - at the same time, increase the power of group 4 by the same amount
    # - apply a small reduction (0.1) in the power for group 2 at t=0  
    # - at the same time, increase the power of group 3 by the same amount
    # - interpolate linearly between t=0 and t=t_end with no correction at t=t_end
    
    base_induction = calc_induction_per_group(group_id, time)
    
    # Calculate interpolation factor (1.0 at t=0, 0.0 at t=t_end)
    if t_end > 0
        interp_factor = max(0.0, (t_end - time) / t_end)
    else
        interp_factor = time == 0 ? 1.0 : 0.0
    end
    
    # Apply corrections based on group
    correction = 0.0
    if group_id == 1
        correction = -0.2 * interp_factor  # Large reduction
    elseif group_id == 4
        correction = +0.2 * interp_factor  # Large increase (balancing group 1)
    elseif group_id == 2
        correction = -0.1 * interp_factor  # Small reduction
    elseif group_id == 3
        correction = +0.1 * interp_factor  # Small increase (balancing group 2)
    end

    rel_power = calc_cp(base_induction) / cp_max + correction
    corrected_induction = calc_induction(rel_power * cp_max)
    return max(0.0, min(1/3, corrected_induction))
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

function calc_induction_matrix(ta, con, time_step, t_end)
    # Create time vector from 0 to t_end with time_step intervals
    time_vector = 0:time_step:t_end
    n_time_steps = length(time_vector)
    n_turbines = size(ta.pos, 1)  # Use ta.pos to get number of turbines
    
    # Initialize matrix: rows = time steps, columns = turbines
    induction_matrix = zeros(Float64, n_time_steps, n_turbines)
    
    # Calculate induction for each turbine at each time step
    for (t_idx, time) in enumerate(time_vector)
        for i in 1:n_turbines
            induction_matrix[t_idx, i] = calc_axial_induction(ta, con, i, time)
        end
    end
    
    return induction_matrix
end

function plot_demand()
    # Create time vector from 0 to t_end with time_step intervals
    time_vector = 0:time_step:t_end

    # Calculate demand for each time point
    demand_values = [calc_demand(t) for t in time_vector]
    
    # Calculate actual induction values for representative turbines from each group
    # (This includes the time-dependent corrections implemented in calc_axial_induction)
    
    # Find representative turbines for each group
    turbine_group1 = 1   # Turbine 1 is in group 1 (gets -0.2 reduction)
    turbine_group2 = 2   # Turbine 2 is in group 2 (gets -0.1 reduction)
    turbine_group3 = 3   # Turbine 3 is in group 3 (gets +0.1 increase)
    turbine_group4 = 7   # Turbine 7 is in group 4 (gets +0.2 increase)
    
    # Calculate induction values using actual calc_axial_induction (with corrections)
    induction_group1 = [calc_axial_induction(ta, con, turbine_group1, t) for t in time_vector]
    induction_group2 = [calc_axial_induction(ta, con, turbine_group2, t) for t in time_vector]
    induction_group3 = [calc_axial_induction(ta, con, turbine_group3, t) for t in time_vector]
    induction_group4 = [calc_axial_induction(ta, con, turbine_group4, t) for t in time_vector]
    
    # Combine all data series
    all_data = [demand_values, induction_group1, induction_group2, induction_group3, induction_group4]
    labels = ["Demand", "Group 1 (-0.2)", "Group 2 (-0.1)", "Group 3 (+0.1)", "Group 4 (+0.2)"]
    
    # Create the plot
    ControlPlots.plot(time_vector, all_data, 
                     xlabel="Time [s]", 
                     ylabel="Demand and Induction [-]", 
                     title="Demand and Induction vs Time for All Groups",
                     labels=labels)
end

function plot_induction_matrix()
    # Calculate induction matrix for all turbines over time
    induction_matrix = calc_induction_matrix(ta, con, time_step, t_end)
    time_vector = 0:time_step:t_end
    n_time_steps = size(induction_matrix, 1)
    
    # Initialize group data arrays
    group_data = [Float64[] for _ in 1:4]  # Arrays for groups 1-4
    
    # Calculate average induction for each group at each time step
    for t_idx in 1:n_time_steps
        group_sums = zeros(4)
        group_counts = zeros(Int, 4)
        
        # Sum induction values for each group
        for turbine in 1:size(ta.pos, 1)
            group_id = FLORIDyn.turbine_group(ta, turbine)
            if 1 <= group_id <= 4
                group_sums[group_id] += induction_matrix[t_idx, turbine]
                group_counts[group_id] += 1
            end
        end
        
        # Calculate averages and store
        for group in 1:4
            if group_counts[group] > 0
                avg_induction = group_sums[group] / group_counts[group]
            else
                avg_induction = 0.0
            end
            push!(group_data[group], avg_induction)
        end
    end
    
    # Create labels with group information
    group_labels = ["Group 1 (-0.2)", "Group 2 (-0.1)", "Group 3 (+0.1)", "Group 4 (+0.2)"]
    
    ControlPlots.plot(time_vector, group_data,
                     xlabel="Time [s]", 
                     ylabel="Axial Induction Factor [-]", 
                     title="Average Axial Induction Factor vs Time by Turbine Group",
                     labels=group_labels)
end

con.induction_data = calc_induction_matrix(ta, con, time_step, t_end)

nothing
