# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

#=
This file contains the main FLORIS wake model execution functions and their specialized helpers.

Functions defined in this file:
- FLORISBuffers(n_pts::Int): Constructor for FLORISBuffers struct
- prepare_rotor_points!: Prepare rotor discretization points with scaling, rotation, and translation
- handle_single_turbine!: Handle the special case when there is only one turbine in the simulation
- setup_computation_buffers!: Initialize and setup computation buffers for multi-turbine wake calculations
- compute_wake_effects!: Compute wake effects for a single upstream turbine on the downstream turbine
- compute_final_wind_shear!: Compute final wind shear effects and effective wind speed for the last turbine
- runFLORIS!: Main orchestrating function that coordinates the FLORIS wake model execution

The runFLORIS! function serves as the main entry point and coordinates execution through the
specialized helper functions. The FLORISBuffers struct is defined in structs_floris.jl.
=#

"""
    FLORISBuffers(n_pts::Int) -> FLORISBuffers

Constructor for FLORISBuffers struct that pre-allocates all necessary arrays for FLORIS computation.

# Arguments
- `n_pts::Int`: Number of rotor discretization points to allocate for

# Returns
- `FLORISBuffers`: Initialized struct with all arrays pre-allocated to size `n_pts`

# Notes
- Result arrays (T_red_arr, T_aTI_arr, T_Ueff, T_weight) are initialized as empty and resized during computation
- This constructor minimizes allocations during FLORIS wake model execution
"""
function FLORISBuffers(n_pts::Int)
    return FLORISBuffers(
        Matrix{Float64}(undef, n_pts, 3),  # tmp_RPs
        Matrix{Float64}(undef, n_pts, 3),  # rotor_pts
        Vector{Float64}(undef, n_pts),     # sig_y
        Vector{Float64}(undef, n_pts),     # sig_z
        Vector{Float64}(undef, n_pts),     # x_0
        Matrix{Float64}(undef, n_pts, 2),  # delta
        Vector{Float64}(undef, n_pts),     # pc_y
        Vector{Float64}(undef, n_pts),     # pc_z
        Vector{Float64}(undef, n_pts),     # cw_y
        Vector{Float64}(undef, n_pts),     # cw_z
        Vector{Float64}(undef, n_pts),     # phi_cw
        Vector{Float64}(undef, n_pts),     # r_cw
        Vector{Bool}(undef, n_pts),        # core
        Vector{Bool}(undef, n_pts),        # nw
        Vector{Bool}(undef, n_pts),        # fw
        Vector{Float64}(undef, n_pts),     # tmp_RPs_r
        Vector{Float64}(undef, n_pts),     # gaussAbs
        Vector{Float64}(undef, n_pts),     # gaussWght
        Vector{Float64}(undef, n_pts),     # exp_y
        Vector{Float64}(undef, n_pts),     # exp_z
        Vector{Bool}(undef, n_pts),        # not_core
        Float64[],                         # T_red_arr (size set per call)
        Float64[],                         # T_aTI_arr (size set per call)
        Float64[],                         # T_Ueff (size 0 or 1)
        Float64[],                         # T_weight (size set per call)
    )
end

"""
    prepare_rotor_points!(buffers::FLORISBuffers, location_t, states_t, d_rotor, floris::Floris)

Prepare rotor discretization points with scaling, rotation, and translation for the last turbine.

This function handles the rotor point discretization, scales them by the rotor diameter,
applies yaw rotation, and translates them to the turbine location. The results are stored
in the buffers to avoid allocations.

# Arguments
- `buffers::FLORISBuffers`: Pre-allocated computation buffers
- `location_t`: Turbine locations matrix
- `states_t`: Turbine states matrix (includes yaw angles)
- `d_rotor`: Rotor diameter array
- `floris::Floris`: FLORIS model parameters

# Returns
- `RPl`: View of transformed rotor points
- `RPw`: Rotor point weights

# Note
This function is **private** and intended for internal use only.
"""
function prepare_rotor_points!(buffers::FLORISBuffers, location_t, states_t, d_rotor, floris::Floris)
    if d_rotor[end] > 0
        RPl, RPw = discretizeRotor(floris.rotor_points)
    else
        RPl = SA[0.0 0.0 0.0]
        RPw = SA[1.0]
    end
    
    # Yaw rotation for last turbine
    tmp_yaw = deg2rad(states_t[end, 2])
    R = SA[cos(tmp_yaw)  sin(tmp_yaw)  0.0;
          -sin(tmp_yaw)  cos(tmp_yaw)  0.0;
           0.0           0.0           1.0]

    # Conservative allocation reduction: scale into rotor_pts buffer (no scaled_RPl alloc),
    # then perform the rotation+translation in-place
    nRP_local = size(RPl, 1)
    # Safety: ensure buffers are large enough before writing to avoid OOB/segfaults
    if size(buffers.rotor_pts, 1) < nRP_local
        error("FLORISBuffers.rotor_pts too small: expected at least $(nRP_local) rows, got $(size(buffers.rotor_pts, 1)).\n" *
              "Ensure create_unified_buffers(.., floris) used the same rotor discretization.")
    end
    @inbounds for i in 1:nRP_local
        buffers.rotor_pts[i, 1] = RPl[i, 1] * d_rotor[end]
        buffers.rotor_pts[i, 2] = RPl[i, 2] * d_rotor[end]
        buffers.rotor_pts[i, 3] = RPl[i, 3] * d_rotor[end]
    end
    scaled_RPl = @view buffers.rotor_pts[1:nRP_local, :]
    
    # In-place rotation + translation to avoid allocations from transpose/mul/broadcast
    @inbounds begin
        # Extract rotation matrix entries (StaticArrays, so this is cheap)
        R11, R12, R13 = R[1,1], R[1,2], R[1,3]
        R21, R22, R23 = R[2,1], R[2,2], R[2,3]
        R31, R32, R33 = R[3,1], R[3,2], R[3,3]
        # Extract translation components once
        tx = location_t[end, 1]
        ty = location_t[end, 2]
        tz = location_t[end, 3]
        for i in 1:nRP_local
            x = scaled_RPl[i, 1]
            y = scaled_RPl[i, 2]
            z = scaled_RPl[i, 3]
            scaled_RPl[i, 1] = R11*x + R12*y + R13*z + tx
            scaled_RPl[i, 2] = R21*x + R22*y + R23*z + ty
            scaled_RPl[i, 3] = R31*x + R32*y + R33*z + tz
        end
    end
    
    return scaled_RPl, RPw
end

"""
    handle_single_turbine!(buffers::FLORISBuffers, RPl, RPw, location_t, set::Settings, windshear, d_rotor)

Handle the special case when there is only one turbine in the simulation.

This function computes the wind shear reduction for a single turbine case and populates
the appropriate buffer arrays. It returns early to avoid the multi-turbine wake calculations.

# Arguments
- `buffers::FLORISBuffers`: Pre-allocated computation buffers
- `RPl`: Rotor discretization points
- `RPw`: Rotor point weights
- `location_t`: Turbine locations matrix
- `set::Settings`: Simulation settings
- `windshear`: Wind shear data
- `d_rotor`: Rotor diameter array

# Returns
- `nothing` (results stored in buffers)

# Note
This function is **private** and intended for internal use only.
"""
function handle_single_turbine!(buffers::FLORISBuffers, RPl, RPw, location_t, set::Settings, 
                               windshear, d_rotor)
    # Avoid allocating RPl[:,3] and the broadcasted division by using a buffer
    nRP_local = size(RPl, 1)
    if length(buffers.tmp_RPs_r) < nRP_local
        error("FLORISBuffers.tmp_RPs_r too small: expected at least $(nRP_local) elements, got $(length(buffers.tmp_RPs_r)).\n" *
              "Ensure create_unified_buffers(.., floris) used the same rotor discretization.")
    end
    if set.shear_mode isa Shear_PowerLaw
        # Power law expects z normalized by hub height; clamp to > 0
        @inbounds for i in 1:nRP_local
            val = RPl[i, 3] / location_t[end, 3]
            buffers.tmp_RPs_r[i] = val > eps() ? val : eps()
        end
    else
        # Interpolation expects absolute height in meters
        @inbounds for i in 1:nRP_local
            buffers.tmp_RPs_r[i] = RPl[i, 3]
        end
    end
    z_view = @view buffers.tmp_RPs_r[1:nRP_local]
    redShear = getWindShearT(set.shear_mode, windshear, z_view)
    # Avoid allocating a view for RPw in the dot product
    acc = 0.0
    @inbounds for i in 1:nRP_local
        acc = muladd(RPw[i], redShear[i], acc)
    end
    T_red_scalar = acc
    # Persist result into buffers as 1-length arrays (optional use by callers)
    resize!(buffers.T_red_arr, 1); buffers.T_red_arr[1] = T_red_scalar
    resize!(buffers.T_aTI_arr, 0)
    resize!(buffers.T_Ueff, 0)
    resize!(buffers.T_weight, 0)
    return nothing
end

"""
    setup_computation_buffers!(buffers::FLORISBuffers, nRP::Int, nT::Int)

Initialize and setup computation buffers for multi-turbine wake calculations.

This function resizes output arrays and creates views of pre-allocated buffers to match
the current rotor discretization size. It ensures all buffers are properly sized before
the main computation loop.

# Arguments
- `buffers::FLORISBuffers`: Pre-allocated computation buffers
- `nRP::Int`: Number of rotor points
- `nT::Int`: Number of turbines

# Returns
- Tuple of buffer views for use in wake calculations

# Note
This function is **private** and intended for internal use only.
"""
function setup_computation_buffers!(buffers::FLORISBuffers, nRP::Int, nT::Int)
    # Initialize outputs in buffers
    resize!(buffers.T_red_arr, nT); fill!(buffers.T_red_arr, 1.0)
    resize!(buffers.T_aTI_arr, max(nT - 1, 0))
    if nT > 1
        fill!(buffers.T_aTI_arr, 0.0)
    end
    resize!(buffers.T_weight, max(nT - 1, 0))
    if nT > 1
        fill!(buffers.T_weight, 0.0)
    end

    # Ensure buffers are properly sized
    if size(buffers.tmp_RPs, 1) < nRP
        error("Buffer tmp_RPs is too small: expected at least $(nRP) rows, got $(size(buffers.tmp_RPs, 1))")
    end
    if length(buffers.cw_y) < nRP
        error("Buffer arrays are too small: expected at least $(nRP) elements, got $(length(buffers.cw_y))")
    end
    
    # Use views of pre-allocated buffers to match the current discretization size exactly
    tmp_RPs = view(buffers.tmp_RPs, 1:nRP, :)
    sig_y = view(buffers.sig_y, 1:nRP)
    sig_z = view(buffers.sig_z, 1:nRP)
    x_0   = view(buffers.x_0, 1:nRP)
    delta = view(buffers.delta, 1:nRP, :)
    pc_y  = view(buffers.pc_y, 1:nRP)
    pc_z  = view(buffers.pc_z, 1:nRP)
    cw_y = view(buffers.cw_y, 1:nRP)
    cw_z = view(buffers.cw_z, 1:nRP)
    phi_cw = view(buffers.phi_cw, 1:nRP)
    r_cw = view(buffers.r_cw, 1:nRP)
    core = view(buffers.core, 1:nRP)
    nw = view(buffers.nw, 1:nRP)
    fw = view(buffers.fw, 1:nRP)
    tmp_RPs_r = view(buffers.tmp_RPs_r, 1:nRP)
    gaussAbs = view(buffers.gaussAbs, 1:nRP)
    gaussWght = view(buffers.gaussWght, 1:nRP)
    exp_y = view(buffers.exp_y, 1:nRP)
    exp_z = view(buffers.exp_z, 1:nRP)
    not_core = view(buffers.not_core, 1:nRP)

    return (tmp_RPs, sig_y, sig_z, x_0, delta, pc_y, pc_z, cw_y, cw_z, phi_cw, r_cw, 
            core, nw, fw, tmp_RPs_r, gaussAbs, gaussWght, exp_y, exp_z, not_core)
end

"""
    compute_wake_effects!(buffers::FLORISBuffers, views, iT::Int, RPl, RPw, location_t, 
                         states_wf, states_t, d_rotor, floris::Floris, nRP::Int)

Compute wake effects for a single upstream turbine on the downstream turbine.

This function calculates the velocity deficit and added turbulence caused by turbine iT
on the last turbine in the array. It handles coordinate transformations, wake variable
calculations, and Gaussian wake modeling.

# Arguments
- `buffers::FLORISBuffers`: Pre-allocated computation buffers
- `views`: Tuple of buffer views from setup_computation_buffers!
- `iT::Int`: Index of the upstream turbine
- `RPl`: Rotor discretization points
- `RPw`: Rotor point weights
- `location_t`: Turbine locations matrix
- `states_wf`: Wind farm states matrix
- `states_t`: Turbine states matrix
- `d_rotor`: Rotor diameter array
- `floris::Floris`: FLORIS model parameters
- `nRP::Int`: Number of rotor points

# Returns
- `nothing`. Results are written into fields of `buffers`:
    - `buffers.T_red_arr`: Velocity reduction factors for each turbine
    - `buffers.T_aTI_arr`: Added turbulence intensity values
    - `buffers.T_Ueff`: Effective wind speeds
    - `buffers.T_weight`: Wake weighting factors

# Note
This function is **private** and intended for internal use only.
"""
function compute_wake_effects!(buffers::FLORISBuffers, views, iT::Int, RPl, RPw, location_t, 
                              states_wf, states_t, d_rotor, floris::Floris, nRP::Int)
    tmp_RPs, sig_y, sig_z, x_0, delta, pc_y, pc_z, cw_y, cw_z, phi_cw, r_cw, 
    core, nw, fw, tmp_RPs_r, gaussAbs, gaussWght, exp_y, exp_z, not_core = views
    
    tmp_phi = size(states_wf,2) == 4 ? angSOWFA2world(states_wf[iT, 4]) :
                                       angSOWFA2world(states_wf[iT, 2])

    # Use pre-allocated array instead of creating new one
    for i in 1:nRP
        for j in 1:3
            tmp_RPs[i, j] = RPl[i, j] - location_t[iT, j]
        end
    end
    
    cos_phi = cos(tmp_phi)
    sin_phi = sin(tmp_phi)
    
    # Apply rotation matrix manually to avoid allocation
    for i in 1:nRP
        x = tmp_RPs[i, 1]
        y = tmp_RPs[i, 2]
        z = tmp_RPs[i, 3]
        tmp_RPs[i, 1] = cos_phi * x + sin_phi * y
        tmp_RPs[i, 2] = -sin_phi * x + cos_phi * y
        tmp_RPs[i, 3] = z
    end

    if tmp_RPs[1, 1] <= 10
        return nothing
    end

    a_val = states_t[iT, 1]
    yaw_deg = states_t[iT, 2]
    yaw = -deg2rad(yaw_deg)
    TI = states_t[iT, 3]
    Ct = calcCt(a_val, yaw_deg)
    TI0 = states_wf[iT, 3]

    # Compute mean_x now, before tmp_RPs is reused for other data
    mean_x = 0.0
    @inbounds for i in 1:nRP
        mean_x += tmp_RPs[i, 1]
    end
    mean_x /= nRP

    # Compute wake variables using in-place API with preallocated buffers
    getVars!(sig_y, sig_z, x_0, delta, pc_y, pc_z, tmp_RPs, Ct, yaw, TI, TI0, floris, d_rotor[iT])
    C_T = Ct

    # Use pre-allocated arrays
    @inbounds for i in 1:nRP
        delta_y = delta[i, 1]
        delta_z = delta[i, 2]
        
        cw_y[i] = tmp_RPs[i, 2] - delta_y
        cw_z[i] = tmp_RPs[i, 3] - delta_z
        phi_cw[i] = atan(cw_z[i], cw_y[i])
        r_cw[i] = sqrt(cw_y[i]^2 + cw_z[i]^2)
        
        pc_y_val = pc_y[i]
        pc_z_val = pc_z[i]
        pc_y_half_cos = 0.5 * pc_y_val * cos(phi_cw[i])
        pc_z_half_sin = 0.5 * pc_z_val * sin(phi_cw[i])
        core[i] = (r_cw[i] < sqrt(pc_y_half_cos^2 + pc_z_half_sin^2)) || (tmp_RPs[i, 1] == 0.0)
        
        nw[i] = tmp_RPs[i, 1] < x_0[i]
        fw[i] = !nw[i]
    end

    # Initialize arrays
    fill!(tmp_RPs_r, 0.0)
    fill!(gaussAbs, 0.0)
    fill!(gaussWght, 1.0)
    
    sqrt_1_minus_CT = sqrt(1 - C_T)
    @inbounds for i in 1:nRP
        if core[i]
            tmp_RPs_r[i] = 1 - sqrt_1_minus_CT
        end
        
        if nw[i]
            gaussAbs[i] = 1 - sqrt_1_minus_CT
        elseif fw[i]
            sig_y_val = sig_y[i]
            sig_z_val = sig_z[i]
            gaussAbs[i] = 1 - sqrt(1 - C_T * cos(yaw) / (8 * sig_y_val * sig_z_val / d_rotor[iT]^2))
        end
    end

    # Handle not_core calculations
    @inbounds for i in 1:nRP
        not_core[i] = !core[i]
    end
    
    any_not_core = false
    @inbounds for i in 1:nRP
        if not_core[i]
            any_not_core = true
            break
        end
    end
    
    if any_not_core
        @inbounds for i in 1:nRP
            if not_core[i]
                sig_y_val = sig_y[i]
                sig_z_val = sig_z[i]
                pc_y_val = pc_y[i]
                pc_z_val = pc_z[i]
                
                y_term = (cw_y[i] - cos(phi_cw[i]) * pc_y_val * 0.5) / sig_y_val
                z_term = (cw_z[i] - sin(phi_cw[i]) * pc_z_val * 0.5) / sig_z_val
                exp_y[i] = exp(-0.5 * y_term^2)
                exp_z[i] = exp(-0.5 * z_term^2)
                gaussWght[i] = exp_y[i] * exp_z[i]
                tmp_RPs_r[i] = gaussAbs[i] * gaussWght[i]
            end
        end
    end

    buffers.T_weight[iT] = sum(gaussWght)
    buffers.T_red_arr[iT] = 1 - dot(RPw, tmp_RPs_r)

    # Added TI
    T_addedTI_tmp = floris.k_fa * (
        a_val^floris.k_fb *
        TI0^floris.k_fc *
        (mean_x / d_rotor[iT])^floris.k_fd
    )

    TIexp = floris.TIexp
    @inbounds for i in 1:nRP
        sig_y_val = sig_y[i]
        sig_z_val = sig_z[i]
        pc_y_val = pc_y[i]
        pc_z_val = pc_z[i]
        
        y_term = (cw_y[i] - cos(phi_cw[i]) * pc_y_val * 0.5) / (TIexp * sig_y_val)
        z_term = (cw_z[i] - sin(phi_cw[i]) * pc_z_val * 0.5) / (TIexp * sig_z_val)
        exp_y[i] = exp(-0.5 * y_term^2)
        exp_z[i] = exp(-0.5 * z_term^2)
    end

    # Avoid allocating a temporary from (exp_y .* exp_z) by manual accumulation
    acc = 0.0
    @inbounds for i in 1:nRP
        acc = muladd(RPw[i], exp_y[i] * exp_z[i], acc)
    end
    buffers.T_aTI_arr[iT] = T_addedTI_tmp * acc
    nothing
end

"""
    compute_final_wind_shear!(buffers::FLORISBuffers, RPl, RPw, location_t, set::Settings, 
                             windshear, tmp_RPs_r, states_wf)

Compute final wind shear effects and effective wind speed for the last turbine.

This function calculates the wind shear reduction for the downstream turbine and computes
the final effective wind speed by combining all wake effects and wind shear.

# Arguments
- `buffers::FLORISBuffers`: Pre-allocated computation buffers
- `RPl`: Rotor discretization points
- `RPw`: Rotor point weights
- `location_t`: Turbine locations matrix
- `set::Settings`: Simulation settings
- `windshear`: Wind shear data
- `tmp_RPs_r`: Temporary buffer for rotor point calculations
- `states_wf`: Wind farm states matrix

# Returns
- `nothing`. Results are written into fields of `buffers`:
    - `buffers.T_red_arr`: Velocity reduction factors for each turbine
    - `buffers.T_aTI_arr`: Added turbulence intensity values
    - `buffers.T_Ueff`: Effective wind speeds
    - `buffers.T_weight`: Wake weighting factors

# Note
This function is **private** and intended for internal use only.
"""
function compute_final_wind_shear!(buffers::FLORISBuffers, RPl, RPw, location_t, set::Settings, 
                                  windshear, tmp_RPs_r, states_wf)
    nRP = size(RPl, 1)
    
    # Compute z for wind shear:
    # - Power law expects normalized height (clamped positive)
    # - Interpolation expects absolute height (meters)
    if set.shear_mode isa Shear_PowerLaw
        @inbounds for i in 1:nRP
            val = RPl[i, 3] / location_t[end, 3]
            tmp_RPs_r[i] = val > eps() ? val : eps()
        end
    else
        @inbounds for i in 1:nRP
            tmp_RPs_r[i] = RPl[i, 3]
        end
    end
    redShear = getWindShearT(set.shear_mode, windshear, tmp_RPs_r)
    buffers.T_red_arr[end] = dot(RPw, redShear)

    T_red = prod(buffers.T_red_arr)
    T_Ueff_scalar = states_wf[end, 1] * T_red
    resize!(buffers.T_Ueff, 1); buffers.T_Ueff[1] = T_Ueff_scalar
    nothing
end

"""
    runFLORIS!(buffers::FLORISBuffers, set::Settings, location_t, states_wf, states_t, d_rotor, 
               floris, windshear::Union{Matrix, WindShear})

Execute the FLORIS (FLOw Redirection and Induction in Steady State) wake model simulation for wind farm analysis.

This is the main orchestrating function that coordinates the FLORIS wake model execution through 
a series of specialized sub-functions. It performs wake analysis using the Gaussian wake model to 
calculate velocity reductions, turbulence intensity additions, and effective wind speeds at turbine 
locations, accounting for wake interactions, rotor discretization, wind shear effects, and turbulence propagation.

# Arguments
- `buffers::FLORISBuffers`: Pre-allocated buffer arrays for computation and output storage (see [`FLORISBuffers`](@ref))
- `set::Settings`: Simulation settings containing configuration options for wind shear modeling
- `location_t`: Matrix of turbine positions [x, y, z] coordinates for each turbine [m]
- `states_wf`: Wind field state matrix containing velocity, direction, and turbulence data
- `states_t`: Turbine state matrix with axial induction factors, yaw angles, and turbulence intensities
- `d_rotor`: Vector of rotor diameters for each turbine [m]
- `floris`: FLORIS model parameters containing wake model coefficients and rotor discretization settings (see [`Floris`](@ref))
- `windshear`: Wind shear profile data for vertical wind speed variation modeling, either a matrix or of type [`WindShear`](@ref)

# Returns
- `nothing`. Results are written into fields of `buffers`:
    - `buffers.T_red_arr`: Velocity reduction factors for each turbine
    - `buffers.T_aTI_arr`: Added turbulence intensity values
    - `buffers.T_Ueff`: Effective wind speeds  
    - `buffers.T_weight`: Wake weighting factors

# Implementation Structure
The function is implemented as a high-level orchestrator that calls specialized sub-functions:

1. **[`prepare_rotor_points!`](@ref)**: Handles rotor discretization, scaling, rotation, and translation
2. **[`handle_single_turbine!`](@ref)**: Optimized path for single turbine simulations (wind shear only)
3. **[`setup_computation_buffers!`](@ref)**: Initializes and configures computation buffers and views
4. **[`compute_wake_effects!`](@ref)**: Core wake computation loop for each upstream turbine
5. **[`compute_final_wind_shear!`](@ref)**: Final wind shear and effective wind speed calculation

# Computational Process
## Single Turbine Case
For single turbine simulations, the function uses an optimized path that only calculates 
wind shear effects without wake interactions, providing significant performance benefits.

## Multi-Turbine Wake Analysis
For multi-turbine wind farms, the function performs:

### Rotor Point Preparation
- Discretizes rotor planes using the isocell algorithm
- Applies yaw rotation transformations and coordinate translations
- Handles both active turbines (d_rotor > 0) and placeholder turbines

### Wake Interaction Calculations
For each upstream turbine affecting downstream turbines:
- **Coordinate Transformations**: Aligns coordinates with wind direction and turbine yaw
- **Wake Variable Calculations**: Computes wake expansion, deflection, and potential core dimensions
- **Velocity Deficit Modeling**: Applies Gaussian wake theory with yaw corrections
- **Turbulence Enhancement**: Calculates added turbulence using empirical correlations
- **Superposition**: Combines effects from multiple upstream turbines

### Final Integration
- Applies vertical wind shear corrections to the downstream turbine
- Combines all wake effects with wind shear for final effective wind speed

# Mathematical Models
The function implements state-of-the-art wake modeling based on:
- **Gaussian Wake Theory**: For velocity deficit calculations with yaw corrections
- **Analytical Deflection Models**: For wake steering in yawed conditions
- **Empirical Turbulence Models**: For wake-added turbulence intensity
- **Linear Superposition**: For combining multiple wake interactions

# Performance Characteristics
- **Computational Complexity**: O(N²) for N turbines due to wake interactions
- **Memory Efficiency**: Uses pre-allocated buffers to avoid runtime allocations
- **Optimization**: Single turbine case bypasses multi-turbine calculations
- **Vectorization**: Leverages SIMD operations where possible

# Notes
- Uses SOWFA (Simulator for Offshore Wind Farm Applications) coordinate conventions
- Supports both research and engineering applications for wind farm optimization
- Requires proper initialization of turbine states and wind field conditions
- Buffer sizes must be compatible with rotor discretization settings

# References
- Bastankhah, M. and Porté-Agel, F. (2016). Experimental and theoretical study of wind turbine wakes in yawed conditions
- Niayifar, A. and Porté-Agel, F. (2016). Analytical modeling of wind farms: A new approach for power prediction

# See Also
- [`prepare_rotor_points!`](@ref): Rotor point preparation and transformation
- [`handle_single_turbine!`](@ref): Single turbine optimization path
- [`setup_computation_buffers!`](@ref): Buffer initialization and setup
- [`compute_wake_effects!`](@ref): Core wake interaction calculations
- [`compute_final_wind_shear!`](@ref): Final wind shear integration
- [`FLORISBuffers`](@ref): Buffer structure documentation
- [`Floris`](@ref): FLORIS model parameters
"""
function runFLORIS!(buffers::FLORISBuffers, set::Settings, location_t, states_wf, states_t, d_rotor, floris, 
                   windshear::Union{Matrix, WindShear})
    # Prepare rotor points (RPl, RPw)
    RPl, RPw = prepare_rotor_points!(buffers, location_t, states_t, d_rotor, floris)
    
    # Handle single turbine case
    if length(d_rotor) == 1
        return handle_single_turbine!(buffers, RPl, RPw, location_t, set, windshear, d_rotor)
    end
    
    # Setup computation buffers
    nRP = size(RPl, 1)
    nT = length(d_rotor)
    views = setup_computation_buffers!(buffers, nRP, nT)
    
    # Main wake computation loop
    for iT in 1:(nT - 1)
        compute_wake_effects!(buffers, views, iT, RPl, RPw, location_t, states_wf, 
                             states_t, d_rotor, floris, nRP)
    end
    
    # Final wind shear computation
    compute_final_wind_shear!(buffers, RPl, RPw, location_t, set, windshear, 
                              views[15], states_wf)  # views[15] is tmp_RPs_r
    
    nothing
end
