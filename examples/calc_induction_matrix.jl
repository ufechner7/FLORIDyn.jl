# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Calculate induction matrix for all turbines over time

using FLORIDyn

const cp_max = 16/27  # Betz limit
const dt = 400

# Forward declarations - these functions are expected to be defined elsewhere
# calc_induction(cp) - calculates induction factor from power coefficient
# t_end - simulation end time (global variable)

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
    initial_demand = 0.4
    final_demand = 0.8
    t1 = 240.0 + dt  # Time to start increasing demand
    t2 = 960.0 + dt  # Time to reach final demand
    if time < t1
        return initial_demand
    elseif time < t2
        return initial_demand + (final_demand - initial_demand) * (time - t1) / (t2 - t1)
    else
        return final_demand
    end
end

"""
    calc_induction_per_group(turbine_group, time)

Calculate base induction factor for a turbine group at a given time.
Assumes no wake effects for simplicity.

# Arguments
- `turbine_group`: Group ID of the turbine
- `time`: Current simulation time [s]

# Returns
- Base induction factor for the group
"""
function calc_induction_per_group(turbine_group, time; scaling = 1.22)
    if USE_MPC
        scaling = 1.247
    end
    # simple example: assume no wakes
    demand = calc_demand(time)
    if USE_MPC
        correction = 1-min(((time - 950)/270)^2, 1)
        demand -= correction*0.016
    end
    induction = calc_induction(demand * scaling * cp_max)
    return induction
end

"""
    calc_axial_induction(ta, con, turbine, time)

Calculate axial induction factor for a specific turbine at a given time.
Includes group-based corrections and time interpolation.

# Arguments
- `ta`: TurbineArray containing turbine positions and configuration
- `con`: Controller configuration object
- `turbine`: Turbine index (1-based)
- `time`: Current simulation time [s]

# Returns
- Axial induction factor for the specified turbine
"""
function calc_axial_induction(ta, con, turbine, time; correction_factor=1.8) # max 1.8
    if ! USE_MPC
        correction_factor = 0.0
    end
    # Check if pre-calculated induction data is available
    if hasfield(typeof(con), :induction_data) && !isnothing(con.induction_data)
        # Use pre-calculated data from con.induction_data
        # First column is time, subsequent columns are turbine data
        time_vector = con.induction_data[:, 1]
        
        # Find the time index that corresponds to the requested time
        time_idx = findfirst(t -> t >= time, time_vector)
        if isnothing(time_idx)
            time_idx = length(time_vector)  # Use last time step if time is beyond range
        end
        
        # Return the pre-calculated induction value for this turbine and time
        # turbine data starts from column 2 (column 1 is time)
        if turbine + 1 <= size(con.induction_data, 2) && time_idx <= size(con.induction_data, 1)
            return con.induction_data[time_idx, turbine + 1]
        end
    end
    
    # Fallback to dynamic calculation if no pre-calculated data is available
    group_id = FLORIDyn.turbine_group(ta, turbine)
    
    # Apply corrections based on turbine group and time with linear interpolation
    # Rules:
    # - apply a large reduction (0.2) in the power for group 1 at t=t1
    # - at the same time, increase the power of group 4 by the same amount
    # - apply a small reduction (0.1) in the power for group 2 at t=t1  
    # - at the same time, increase the power of group 3 by the same amount
    # - interpolate linearly between t=0 and t=t_end with no correction at t=t2
    
    base_induction = calc_induction_per_group(group_id, time)
    t1 = 240.0 + dt  # Time to start increasing demand
    t2 = 960.0 + dt  # Time to reach final demand

    # Calculate interpolation factor
    # 1.0 at t=t1 (full correction), 0.0 at t=t2 (no correction)
    # Linear interpolation between t1 and t2
    if time <= t1
        interp_factor = 1.0  # Full correction before and at t1
    elseif time >= t2
        interp_factor = 0.0  # No correction at and after t2
    else
        # Linear interpolation between t1 and t2
        interp_factor = ((t2 - time) / (t2 - t1))
    end
    
    # Apply corrections based on group
    correction = 0.0
    if group_id == 1
        correction = -0.13 * interp_factor  # Large reduction
    elseif group_id == 4
        correction = +0.2 * interp_factor  # Large increase (balancing group 1)
    elseif group_id == 2
        correction = -0.1 * interp_factor  # Small reduction
    elseif group_id == 3
        correction = -0.00 * interp_factor  # Small increase (balancing group 2)
    end
    correction *= correction_factor  # Apply overall correction factor if needed

    rel_power = calc_cp(base_induction) / cp_max + correction
    corrected_induction = calc_induction(rel_power * cp_max)
    return max(0.0, min(1/3, corrected_induction))
end

"""
    calc_induction_matrix(ta, con, time_step, t_end)

Calculate a matrix of axial induction factors for all turbines over time.

# Arguments
- `ta`: TurbineArray containing turbine positions and configuration
- `con`: Controller configuration object
- `time_step`: Time step for the simulation [s]
- `t_end`: End time of the simulation [s]

# Returns
- `induction_matrix`: Matrix where:
  - First column contains time values
  - Subsequent columns contain induction factors for each turbine
  - Rows represent time steps
"""
function calc_induction_matrix(ta, con, time_step, t_end)
    # Create time vector from 0 to t_end with time_step intervals
    time_vector = 0:time_step:t_end
    n_time_steps = length(time_vector)
    n_turbines = size(ta.pos, 1)  # Use ta.pos to get number of turbines
    
    # Initialize matrix: rows = time steps, columns = time + turbines
    # First column is time, subsequent columns are turbine induction factors
    induction_matrix = zeros(Float64, n_time_steps, n_turbines + 1)
    
    # Fill the first column with time values
    induction_matrix[:, 1] = collect(time_vector)
    
    # Calculate induction for each turbine at each time step (columns 2 onwards)
    for (t_idx, time) in enumerate(time_vector)
        for i in 1:n_turbines
            induction_matrix[t_idx, i + 1] = calc_axial_induction(ta, con, i, time)
        end
    end
    
    return induction_matrix
end
