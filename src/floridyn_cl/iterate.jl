# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    iterateOPs!(::IterateOPs_basic, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn)

Advance operational points through the wind field using time-marching dynamics.

This function implements the basic time-stepping algorithm for operational points (OPs) in the 
FLORIDyn simulation. It advances all operational points through space based on wind velocities 
and wake deflection effects, maintaining proper ordering and updating turbine states accordingly.

# Arguments
- `::IterateOPs_basic`: Dispatch type indicating the basic iteration algorithm
- `wf::WindFarm`: Wind farm object containing all turbine and operational point states. See [`WindFarm`](@ref)
  - `wf.States_OP`: Matrix of operational point states [x, y, z, `dw_pos`, `cw_x`, `cw_z`, ...]
  - `wf.States_T`: Matrix of turbine states 
  - `wf.States_WF`: Matrix of wind field states [velocity, direction, `turbulence_intensity`, ...]
  - `wf.StartI`: Starting indices for each turbine's operational points
  - `wf.nT`: Number of turbines
  - `wf.nOP`: Number of operational points per turbine
  - `wf.D`: Rotor diameters [m]
- `sim::Sim`: Simulation configuration object. See [`Sim`](@ref)
  - `sim.time_step`: Time step size [s]
  - `sim.dyn.advection`: Advection scaling factor
- `floris::Floris`: FLORIS model parameters for wake deflection calculations. See [`Floris`](@ref)
- `floridyn::FloriDyn`: FLORIDyn model parameters. See [`FloriDyn`](@ref)

# Returns
- nothing

# Algorithm
The function performs the following steps for each time iteration:

## 1. State Preservation
Saves the current turbine operational point states before advancement to preserve boundary conditions.

## 2. Spatial Advancement
### Downwind Movement
Advances operational points downstream based on local wind velocity:
```
step_dw = Δt × U × advection_factor
OP_dw_position += step_dw
```

### Crosswind Deflection
Calculates wake centerline deflection using FLORIS model and updates crosswind positions:
```
deflection = centerline(States_OP, States_T, States_WF, floris, D)
OP_cw_position = deflection
```

## 3. Coordinate System Transformation
Converts local wind-aligned movements to world coordinates using wind direction:
```
φ = angSOWFA2world(wind_direction)
x_world += cos(φ) × step_dw - sin(φ) × step_cw_x
y_world += sin(φ) × step_dw + cos(φ) × step_cw_x
z_world += step_cw_z
```

## 4. Temporal Shifting
Uses circular shifting to advance the time history:
- Shifts all state matrices by one time step
- Initializes new operational points with saved turbine states
- Maintains temporal continuity of the simulation

## 5. Spatial Ordering
Ensures operational points remain ordered by downstream position for each turbine:
- Sorts operational points by downstream position (`States_OP[:, 4]`)
- Maintains consistency across all state matrices

# Mathematical Description
The coordinate transformation from wind-aligned to world coordinates follows:
```
[x']   [cos(φ)  -sin(φ)] [step_dw ]
[y'] = [sin(φ)   cos(φ)] [step_cw_x]
```
where `φ` is the wind direction angle in world coordinates.

# Notes
- The function modifies the wind farm object in-place (indicated by the `!` suffix)
- Temporal shifting maintains a moving window of operational point history
- Spatial ordering ensures downstream distance monotonicity for wake calculations
- The algorithm handles both advection and deflection physics simultaneously
- Coordinate transformations account for SOWFA wind direction conventions via [`angSOWFA2world`](@ref)
"""
@views function iterateOPs!(::IterateOPs_basic, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn)
    # Save turbine OPs
    tmpOPStates = copy(wf.States_OP[wf.StartI, :])
    tmpTStates  = copy(wf.States_T[wf.StartI, :])
    tmpWFStates = copy(wf.States_WF[wf.StartI, :])

    # Shift states
    # Downwind step
    step_dw = sim.time_step .* sim.dyn.advection .* wf.States_WF[:, 1] 
    wf.States_OP[:, 4] .+= step_dw

    # Crosswind step
    deflection = centerline(wf.States_OP, wf.States_T, wf.States_WF, floris, wf.D[1])
    step_cw = deflection .- wf.States_OP[:, 5:6]
    wf.States_OP[:, 5:6] .= deflection

    # World coordinate system adjustment
    @inbounds @simd for i in axes(wf.States_OP,1)
        phiwi = angSOWFA2world(wf.States_WF[i, 2])  # Convert wind direction to world coordinates
        cphiwi = cos(phiwi)
        sphiwi = sin(phiwi)
        ai = step_cw[i, 1]
        sdwi = step_dw[i]
        wf.States_OP[i, 1] += cphiwi * sdwi - sphiwi * ai
        wf.States_OP[i, 2] += sphiwi * sdwi + cphiwi * ai
        wf.States_OP[i, 3] += step_cw[i, 2]
    end

    # Circshift & init first OPs
    # OPs
    wf.States_OP = circshift(wf.States_OP, (1, 0))
    wf.States_OP[wf.StartI, :] .= tmpOPStates

    # Turbines
    wf.States_T = circshift(wf.States_T, (1, 0))
    wf.States_T[wf.StartI, :] .= tmpTStates

    # Wind Farm
    wf.States_WF = circshift(wf.States_WF, (1, 0))
    wf.States_WF[wf.StartI, :] .= tmpWFStates

    # Check if OPs are in order
    buf = zeros(Int, size(wf.States_OP, 1))
    for iT in 1:wf.nT
        inds = wf.StartI[iT]:(wf.StartI[iT] + wf.nOP - 1)
        indOP = buf[inds]
        sortperm!(indOP,wf.States_OP[inds, 4])
        if ! issorted(indOP)  # check if already sorted
            wf.States_OP[inds, :] .= wf.States_OP[inds[indOP], :]
            wf.States_T[inds, :]  .= wf.States_T[inds[indOP], :]
            wf.States_WF[inds, :] .= wf.States_WF[inds[indOP], :]
        end
    end
    return nothing
end

"""
    IterateOPsBuffers

A struct containing pre-allocated buffers for allocation-free execution of iterateOPs!.

This struct eliminates all allocations during the operational point iteration by 
pre-allocating all necessary temporary arrays. It should be created once and reused
across multiple calls to iterateOPs! for maximum performance.

# Fields
- `tmpOPStates::Matrix{Float64}`: Buffer for saving turbine operational point states  
- `tmpTStates::Matrix{Float64}`: Buffer for saving turbine states
- `tmpWFStates::Matrix{Float64}`: Buffer for saving wind farm states
- `step_dw::Vector{Float64}`: Buffer for downwind step calculations
- `deflection::Matrix{Float64}`: Buffer for centerline deflection calculations  
- `step_cw::Matrix{Float64}`: Buffer for crosswind step calculations
- `temp_states_op::Matrix{Float64}`: Temporary buffer for States_OP circular shifting
- `temp_states_t::Matrix{Float64}`: Temporary buffer for States_T circular shifting  
- `temp_states_wf::Matrix{Float64}`: Temporary buffer for States_WF circular shifting
- `sort_buffer::Vector{Int}`: Buffer for sorting operational points

# Constructor
    IterateOPsBuffers(wf::WindFarm)

Creates buffers appropriately sized for the given WindFarm object.
"""
mutable struct IterateOPsBuffers
    tmpOPStates::Matrix{Float64}     # For saving initial turbine OP states (nT × states_op)
    tmpTStates::Matrix{Float64}      # For saving initial turbine states (nT × states_t)  
    tmpWFStates::Matrix{Float64}     # For saving initial wind farm states (nT × states_wf)
    step_dw::Vector{Float64}         # For downwind step calculation (total_ops)
    deflection::Matrix{Float64}      # For centerline deflection (total_ops × 2)
    step_cw::Matrix{Float64}         # For crosswind step (total_ops × 2)
    temp_states_op::Matrix{Float64}  # Temporary for circular shift (total_ops × states_op)
    temp_states_t::Matrix{Float64}   # Temporary for circular shift (total_ops × states_t) 
    temp_states_wf::Matrix{Float64}  # Temporary for circular shift (total_ops × states_wf)
    sort_buffer::Vector{Int}         # For sorting permutations (total_ops)
end

function IterateOPsBuffers(wf::WindFarm)
    n_turbines = size(wf.StartI, 2)
    total_ops = size(wf.States_OP, 1)
    n_states_op = size(wf.States_OP, 2)
    n_states_t = size(wf.States_T, 2)
    n_states_wf = size(wf.States_WF, 2)
    
    return IterateOPsBuffers(
        Matrix{Float64}(undef, n_turbines, n_states_op),    # tmpOPStates
        Matrix{Float64}(undef, n_turbines, n_states_t),     # tmpTStates
        Matrix{Float64}(undef, n_turbines, n_states_wf),    # tmpWFStates
        Vector{Float64}(undef, total_ops),                  # step_dw
        Matrix{Float64}(undef, total_ops, 2),               # deflection
        Matrix{Float64}(undef, total_ops, 2),               # step_cw
        Matrix{Float64}(undef, total_ops, n_states_op),     # temp_states_op
        Matrix{Float64}(undef, total_ops, n_states_t),      # temp_states_t
        Matrix{Float64}(undef, total_ops, n_states_wf),     # temp_states_wf
        Vector{Int}(undef, total_ops)                       # sort_buffer
    )
end

"""
    iterateOPs!(::IterateOPs_basic, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn, buffers::IterateOPsBuffers)

Allocation-free version of iterateOPs! that uses pre-allocated buffers.

This is the high-performance, zero-allocation version of the operational point iteration algorithm.
All temporary arrays are pre-allocated in the buffers parameter to eliminate runtime allocations.

# Arguments  
- `::IterateOPs_basic`: Dispatch type indicating the basic iteration algorithm
- `wf::WindFarm`: Wind farm object (same as standard version)
- `sim::Sim`: Simulation configuration object (same as standard version)  
- `floris::Floris`: FLORIS model parameters (same as standard version)
- `floridyn::FloriDyn`: FLORIDyn model parameters (same as standard version)
- `buffers::IterateOPsBuffers`: Pre-allocated buffers for all temporary calculations

# Returns
- nothing

# Performance Notes
- Zero allocations during execution (after initial buffer setup)
- Suitable for performance-critical applications and benchmarking
- Buffers can be reused across multiple calls for maximum efficiency
"""
@views function iterateOPs!(::IterateOPs_basic, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn, buffers::IterateOPsBuffers)
    # Save turbine OPs using pre-allocated buffers
    @inbounds for i in 1:size(wf.StartI, 2)
        start_idx = wf.StartI[1, i]
        @views buffers.tmpOPStates[i, :] .= wf.States_OP[start_idx, :]
        @views buffers.tmpTStates[i, :] .= wf.States_T[start_idx, :]
        @views buffers.tmpWFStates[i, :] .= wf.States_WF[start_idx, :]
    end

    # Compute downwind step without allocation
    @inbounds @simd for i in 1:size(wf.States_OP, 1)
        buffers.step_dw[i] = sim.time_step * sim.dyn.advection * wf.States_WF[i, 1]
    end
    
    # Apply downwind step
    @inbounds @simd for i in 1:size(wf.States_OP, 1)
        wf.States_OP[i, 4] += buffers.step_dw[i]
    end

    # Crosswind step - compute deflection using in-place function
    centerline!(buffers.deflection, wf.States_OP, wf.States_T, wf.States_WF, floris, wf.D[1])
    
    # Compute crosswind step without allocation
    @inbounds for i in 1:size(wf.States_OP, 1)
        buffers.step_cw[i, 1] = buffers.deflection[i, 1] - wf.States_OP[i, 5]
        buffers.step_cw[i, 2] = buffers.deflection[i, 2] - wf.States_OP[i, 6]
    end
    
    # Update crosswind positions
    @inbounds for i in 1:size(wf.States_OP, 1)
        wf.States_OP[i, 5] = buffers.deflection[i, 1]
        wf.States_OP[i, 6] = buffers.deflection[i, 2]
    end

    # World coordinate system adjustment
    @inbounds @simd for i in axes(wf.States_OP,1)
        phiwi = angSOWFA2world(wf.States_WF[i, 2])
        cphiwi = cos(phiwi)
        sphiwi = sin(phiwi)
        ai = buffers.step_cw[i, 1]
        sdwi = buffers.step_dw[i]
        wf.States_OP[i, 1] += cphiwi * sdwi - sphiwi * ai
        wf.States_OP[i, 2] += sphiwi * sdwi + cphiwi * ai
        wf.States_OP[i, 3] += buffers.step_cw[i, 2]
    end

    # Manual circular shift for States_OP
    _circshift_and_restore!(wf.States_OP, buffers.tmpOPStates, wf.StartI, buffers.temp_states_op)

    # Manual circular shift for States_T
    _circshift_and_restore!(wf.States_T, buffers.tmpTStates, wf.StartI, buffers.temp_states_t)

    # Manual circular shift for States_WF
    _circshift_and_restore!(wf.States_WF, buffers.tmpWFStates, wf.StartI, buffers.temp_states_wf)

    # Check if OPs are in order and reorder if necessary
    _reorder_ops!(wf, buffers)
    
    return nothing
end

"""
    _circshift_and_restore!(data, initial_states, start_indices, buffer)

Perform an in-place circular shift down by one row and restore initial states.
"""
function _circshift_and_restore!(data::AbstractMatrix, initial_states::AbstractMatrix, start_indices::AbstractMatrix, buffer::AbstractMatrix)
    @inbounds for j in 1:size(data, 2)
        # Copy to buffer
        for i in 1:size(data, 1)
            buffer[i, j] = data[i, j]
        end
        
        # Shift down by one row
        data[1, j] = buffer[end, j]
        for i in 2:size(data, 1)
            data[i, j] = buffer[i-1, j]
        end
    end
    
    # Restore initial states
    @inbounds for i in 1:size(start_indices, 2)
        start_idx = start_indices[1, i]
        @views data[start_idx, :] .= initial_states[i, :]
    end
end

"""
    _reorder_ops!(wf, buffers)

Check if operational points are in order and reorder them if not, using pre-allocated buffers.
"""
function _reorder_ops!(wf::WindFarm, buffers::IterateOPsBuffers)
    @inbounds for iT in 1:wf.nT
        start_idx = wf.StartI[1, iT]
        inds = start_idx:(start_idx + wf.nOP - 1)
        
        # Use a view of the sort buffer for this turbine's range
        sort_perm_view = view(buffers.sort_buffer, 1:wf.nOP)
        
        # Get downstream positions and sort
        downstream_positions = view(wf.States_OP, inds, 4)
        sortperm!(sort_perm_view, downstream_positions)
        
        # Check if reordering is needed
        if !issorted(sort_perm_view)
            # Use views of buffers for reordering
            op_view = view(wf.States_OP, inds, :)
            t_view = view(wf.States_T, inds, :)
            wf_view = view(wf.States_WF, inds, :)
            
            temp_op_view = view(buffers.temp_states_op, 1:wf.nOP, :)
            temp_t_view = view(buffers.temp_states_t, 1:wf.nOP, :)
            temp_wf_view = view(buffers.temp_states_wf, 1:wf.nOP, :)

            # Save current state to buffers before reordering
            temp_op_view .= op_view
            temp_t_view .= t_view
            temp_wf_view .= wf_view
            
            # Reorder based on the sort permutation
            @inbounds for i in 1:wf.nOP
                op_view[i, :] .= view(temp_op_view, sort_perm_view[i], :)
                t_view[i, :] .= view(temp_t_view, sort_perm_view[i], :)
                wf_view[i, :] .= view(temp_wf_view, sort_perm_view[i], :)
            end
        end
    end
end
