# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Calculate induction matrix for all turbines over time

# This file is included by calc_induction.jl and mpc.jl
# It provides the functions:
# - calc_cp(induction)
# - calc_induction(cp)
# - calc_demand(vis, time)
# - calc_induction_per_group(vis, turbine_group, time)
# - calc_induction_matrix(vis, ta, time_step, t_end)
# - calc_axial_induction(vis, ta, turbine, time)

using FLORIDyn

const cp_max = 16/27  # Betz limit
const BETZ_INDUCTION = 1/3

function calc_cp(induction)
    return 4 * induction * (1 - induction)^2
end

"""
    calc_induction(cp)

Calculates the axial induction factor from the power coefficient.

Since we know the induction is always in [0, BETZ_INDUCTION], we can use
a simplified approach that only searches in the physical range.

# Arguments
- `cp::Float64`: Power coefficient (typically in range [0, 16/27])

# Returns
- `Float64`: Axial induction factor in range [0, 1/3]

# Notes
The function solves cp = 4*a*(1-a)^2 for the induction factor 'a'.
"""
function calc_induction(cp)
    # Handle edge cases
    if cp <= 0
        return 0.0
    elseif cp >= cp_max
        return BETZ_INDUCTION
    end
    
    # Since we know a ∈ [0, BETZ_INDUCTION], use simplified bisection
    # Define the function f(a) = 4*a*(1-a)² - cp
    f(a) = 4 * a * (1 - a)^2 - cp
    
    # Only search in the physical range [0, BETZ_INDUCTION]
    a_low = 0.0
    a_high = BETZ_INDUCTION
    
    # Simplified bisection method
    tolerance = 1e-12
    max_iterations = 50  # Reduced since smaller search space
    
    for i in 1:max_iterations
        a_mid = (a_low + a_high) / 2
        f_mid = f(a_mid)
        
        if abs(f_mid) < tolerance || (a_high - a_low) < tolerance
            return a_mid
        end
        
        # Since f is monotonic in [0, BETZ_INDUCTION], use simple comparison
        if f_mid > 0
            a_high = a_mid
        else
            a_low = a_mid
        end
    end
    
    return (a_low + a_high) / 2
end

function calc_demand(vis::Vis, time)
    if USE_STEP
        # Example: step demand profile
        if time < 200 + vis.t_skip
            return 0.001
        else
            return 0.999
        end
    else
        initial_demand = 0.4
        final_demand = 0.8
        t1 = vis.t_skip + T_START  # Time to start increasing demand
        t2 = vis.t_skip + T_END    # Time to reach final demand
        if time < t1
            return initial_demand
        elseif time < t2
            return initial_demand + (final_demand - initial_demand) * (time - t1) / (t2 - t1)
        else
            return final_demand
        end
    end
end

function calc_wind(vis::Vis, time)
    local wind
    low_wind = 6.0
    high_wind = 8.2
    t1 = vis.t_skip + T_START  # Time to start increased wind speed
    t2 = vis.t_skip + T_END    # Time to stop  increased wind speed
    if time < t1
        wind = low_wind
    elseif time < t2
        wind = high_wind
    else
        wind = low_wind
    end

    return wind
end

function calc_vel(vis::Vis, start_time::Real, end_time::Real)
    local vel
    dt = 4
    low_wind = 6.0
    high_wind = 8.2
    # Use absolute times based on start_time
    t1 = start_time + vis.t_skip + T_START  # Absolute time to start increased wind speed
    t2 = start_time + vis.t_skip + T_END    # Absolute time to stop increased wind speed
    # Ensure the velocity table extends to the end of simulation
    vel = [start_time low_wind
            t1-dt low_wind
            t1 high_wind
            t2-dt high_wind
            t2 low_wind
            end_time low_wind]
    return vel
end

"""
    calc_induction_per_group(turbine_group, time; scaling = 1.22)

Calculate base induction factor for a turbine group at a given time.

This function computes the baseline axial induction factor for a turbine group
based on time-varying power demand. The calculation assumes no wake effects
and provides the foundation for group-based wind farm control.

# Arguments
- `turbine_group`: Group ID of the turbine (typically 1-4)
- `time`: Current simulation time [s]
- `scaling=1.22`: Power demand scaling factor (1.247 if `USE_TGC` is enabled)

# Returns
- `Float64`: Base axial induction factor in range [0, 1/3]

# Implementation Details
The function:
1. Calculates time-varying power demand using [`calc_demand`](@ref)
2. Applies TGC-specific corrections if `USE_TGC` is enabled
3. Converts demand to induction factor via power coefficient relationship
4. Assumes uniform behavior across all turbines in the same group

# Notes
This baseline induction can be further modified by [`calc_axial_induction`](@ref)
to include group-specific corrections and wake interactions.

# See Also
- [`calc_demand`](@ref): Time-varying power demand calculation
- [`calc_induction`](@ref): Power coefficient to induction conversion
- [`calc_axial_induction`](@ref): Turbine-specific induction with corrections
"""
function calc_induction_per_group(vis::Vis, turbine_group, time; scaling = 1.22)
    if USE_TGC
        scaling = 1.247
    elseif USE_STEP
        scaling = 1.0
    else
    end
    # simple example: assume no wakes
    demand = calc_demand(vis, time)
    if USE_TGC
        correction = 1-min(((time - 950)/270)^2, 1)
        demand -= correction*0.016
    end
    induction = calc_induction(demand * scaling * cp_max)
    return induction
end

"""
    calc_axial_induction(vis::Vis, ta, turbine, time; correction_factor=1.8)

Calculate axial induction factor for a specific turbine at a given time.
Includes group-based corrections and time interpolation.

This function computes the turbine-specific axial induction factor by applying
time-interpolated corrections to the baseline induction from [`calc_induction_per_group`](@ref).
The corrections implement a group-based control strategy that redistributes power
demand across turbine groups during the transition period between T_START and T_END.

# Arguments
- `vis::Vis`: [`Vis`](@ref) object containing simulation settings (e.g., `t_skip`)
- `ta`: TurbineArray containing turbine positions and configuration
- `turbine`: Turbine index (1-based)
- `time`: Current simulation time [s]
- `correction_factor=1.8`: Maximum correction scaling factor (set to 0.0 if `USE_TGC` is disabled)

# Returns
- `Float64`: Axial induction factor for the specified turbine, clamped to [0, 1/3]

# Implementation Details
The function applies time-interpolated power corrections based on turbine group:
- **Group 1**: -0.13 correction (large reduction)
- **Group 2**: -0.10 correction (small reduction)
- **Group 3**: -0.00 correction (no change)
- **Group 4**: +0.20 correction (large increase, balancing group 1)

Interpolation schedule:
- `time ≤ T_START`: Full correction applied (interpolation factor = 1.0)
- `T_START < time < T_END`: Linear interpolation from full to no correction
- `time ≥ T_END`: No correction applied (interpolation factor = 0.0)

The corrections are scaled by `correction_factor` and converted to induction factors
while respecting the Betz limit (maximum induction = 1/3).

# Global Constants Used
- `T_START`: Time offset to start ramping corrections (relative to `vis.t_skip`)
- `T_END`: Time offset to reach zero correction (relative to `vis.t_skip`)
- `USE_TGC`: Flag to enable/disable turbine group control corrections

# See Also
- [`calc_induction_per_group`](@ref): Baseline induction calculation per group
- [`calc_induction`](@ref): Power coefficient to induction conversion
- [`calc_cp`](@ref): Induction to power coefficient conversion
"""
function calc_axial_induction(vis::Vis, ta, turbine, time; correction_factor=1.8) # max 1.8
    if ! USE_TGC
        correction_factor = 0.0
    end
    group_id = FLORIDyn.turbine_group(ta, turbine)
    
    # Apply corrections based on turbine group and time with linear interpolation
    # Rules:
    # - apply a large reduction (0.2) in the power for group 1 at t=t1
    # - at the same time, increase the power of group 4 by the same amount
    # - apply a small reduction (0.1) in the power for group 2 at t=t1  
    # - at the same time, increase the power of group 3 by the same amount
    # - interpolate linearly between t=0 and t=t_end with no correction at t=t2
    
    base_induction = calc_induction_per_group(vis, group_id, time)
    t1 = vis.t_skip + T_START  # Time to start increasing demand
    t2 = vis.t_skip + T_END    # Time to reach final demand

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
    return max(0.0, min(BETZ_INDUCTION, corrected_induction))
end

"""
    calc_induction_matrix(vis::Vis, ta, time_step, t_end)

Calculate a matrix of axial induction factors for all turbines over time.

# Arguments
- `vis::Vis`: [`Vis`](@ref) containing simulation settings (e.g., `t_skip`)
- `ta`: TurbineArray containing turbine positions and configuration
- `time_step`: Time step for the simulation [s]
- `t_end`: End time of the simulation [s]

# Returns
- `induction_matrix`: Matrix where:
  - First column contains time values
  - Subsequent columns contain induction factors for each turbine
  - Rows represent time steps
"""
function calc_induction_matrix(vis::Vis, ta, time_step, t_end)
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
            induction_matrix[t_idx, i + 1] = calc_axial_induction(vis, ta, i, time)
        end
    end
    
    return induction_matrix
end
