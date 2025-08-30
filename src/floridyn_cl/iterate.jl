# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

#= iterate.jl - Observation point iteration for wind farm simulations

This file contains functions for advancing observation points through the wind field
in wind farm wake modeling simulations. The observation points track wake evolution
through space and time.

Functions and structs defined in this file:
- IterateOPsBuffers: Struct with pre-allocated buffers for allocation-free execution
- IterateOPsBuffers(wf): Constructor for buffer struct
- iterateOPs!: Generic function for observation point iteration (multiple dispatch)
- iterateOPs!(::IterateOPs_basic, wf, sim, floris, floridyn, buffers): High-performance basic iteration algorithm
- _circshift_and_restore!(data, initial_states, start_indices, buffer): In-place circular shift with state restoration
- _reorder_ops!(wf, buffers): Reorder observation points to maintain downstream ordering

The main functionality handles:
- Downwind advection of observation points based on local wind velocity
- _reorder_ops!(wf, buffers): Reorder operational points to maintain downstream ordering

The main functionality handles:
- Downwind advection of operational points based on local wind velocity
>>>>>>> f37a98a0 (Fix plots and improve code organization)
- Crosswind deflection due to wake-induced effects
- Coordinate transformation to world coordinates
- Temporal advancement through circular shifting
- Spatial reordering to maintain proper downstream sequencing =#

"""
    IterateOPsBuffers

A struct containing pre-allocated buffers for allocation-free execution of iterateOPs!.

This struct eliminates all allocations during the observation point iteration by 
pre-allocating all necessary temporary arrays. It should be created once and reused
across multiple calls to iterateOPs! for maximum performance.

# Fields
- `tmpOPStates::Matrix{Float64}`: Buffer for saving turbine observation point states  
- `tmpTStates::Matrix{Float64}`: Buffer for saving turbine states
- `tmpWFStates::Matrix{Float64}`: Buffer for saving wind farm states
- `step_dw::Vector{Float64}`: Buffer for downwind step calculations
- `deflection::Matrix{Float64}`: Buffer for centerline deflection calculations  
- `step_cw::Matrix{Float64}`: Buffer for crosswind step calculations
- `temp_states_op::Matrix{Float64}`: Temporary buffer for States_OP circular shifting
- `temp_states_t::Matrix{Float64}`: Temporary buffer for States_T circular shifting  
- `temp_states_wf::Matrix{Float64}`: Temporary buffer for States_WF circular shifting
- `sort_buffer::Vector{Int}`: Buffer for sorting observation points

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
    iterateOPs!(iterate_mode::IterateOPs_model, wf::WindFarm, sim::Sim, floris::Floris, 
                floridyn::FloriDyn, buffers::IterateOPsBuffers) -> Nothing

Advance observation points through the wind field using the specified iteration strategy.

This function family implements different algorithms for moving observation points (OPs) 
through space and time, which is essential for accurate wake propagation modeling in 
wind farm simulations. The choice of iteration method affects computational efficiency, 
numerical stability, and physical accuracy.

# Summary
The function modifies the following WindFarm fields:
- `wf.States_OP`: Updates observation point positions and states through temporal advancement
- `wf.States_T`: Updates turbine states through circular shifting and temporal evolution  
- `wf.States_WF`: Updates wind field states through circular shifting and temporal evolution

# Input/ Output Arguments
- `wf::WindFarm`: Wind farm object containing turbine and observation point data

# Input Arguments
- `iterate_mode::IterateOPs_model`: Iteration strategy (e.g., [`IterateOPs_basic`](@ref), [`IterateOPs_average`](@ref))
- `sim::Sim`: Simulation configuration with time-stepping parameters
- `floris::Floris`: FLORIS model parameters for wake calculations
- `floridyn::FloriDyn`: FLORIDyn model parameters for wake dynamics
- `buffers::IterateOPsBuffers`: Pre-allocated buffers for allocation-free execution

# Algorithm Overview
1. **State Preservation**: Save initial turbine observation point states
2. **Downwind Advection**: Move OPs downstream based on local wind velocity
3. **Crosswind Deflection**: Apply wake-induced lateral deflection using centerline calculations
4. **Coordinate Transformation**: Convert to world coordinates using wind direction
5. **Temporal Advancement**: Perform circular shifting to advance time steps
6. **Spatial Reordering**: Maintain downstream position ordering of observation points

# Available Methods
- `iterateOPs!(::IterateOPs_basic, ...)`: Basic time-stepping with simple advection

# Notes
- Different iteration strategies provide trade-offs between accuracy and computational cost
"""
function iterateOPs! end

"""
    iterateOPs!(::IterateOPs_basic, wf::WindFarm, sim::Sim, floris::Floris, 
                floridyn::FloriDyn, buffers::IterateOPsBuffers) -> Nothing

This is the high-performance version of the observation point iteration algorithm.

# Buffer Arguments
- `buffers::IterateOPsBuffers`: Pre-allocated buffers for all temporary calculations. See [`IterateOPsBuffers`](@ref).

# Input/ Output Arguments
- `wf::WindFarm`: Wind farm object (same as standard version)

# Input Arguments  
- `::IterateOPs_basic`: Dispatch type indicating the basic iteration algorithm
- `sim::Sim`: Simulation configuration object (same as standard version)  
- `floris::Floris`: FLORIS model parameters (same as standard version)
- `floridyn::FloriDyn`: FLORIDyn model parameters (same as standard version)

# Returns
- nothing
"""
@views function iterateOPs!(::IterateOPs_basic, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn, 
                            buffers::IterateOPsBuffers)
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
    iterateOPs!(::IterateOPs_average, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn, buffers::IterateOPsBuffers) -> Nothing

Observation point iteration with wind field state averaging.

This method mirrors the legacy MATLAB `iterateOPs` implementation where the wind
field (`States_WF`) is advanced using an exponential moving average of the previous
and newly shifted states controlled by weights `sim.dyn.op_iter_weights`.

Weights convention:
- `wNew = sim.dyn.op_iter_weights[1]` (weight for newly shifted state)
- `wOld = sim.dyn.op_iter_weights[2]` (weight for existing state)

The rest of the algorithm (advection, crosswind deflection, coordinate transform,
reordering) matches the basic iteration behavior.
"""
@views function iterateOPs!(::IterateOPs_average, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn, 
                             buffers::IterateOPsBuffers)
    # Save turbine OP, turbine, and wind-field states at start indices
    @inbounds for i in 1:size(wf.StartI, 2)
        start_idx = wf.StartI[1, i]
        @views buffers.tmpOPStates[i, :] .= wf.States_OP[start_idx, :]
        @views buffers.tmpTStates[i, :]   .= wf.States_T[start_idx, :]
        @views buffers.tmpWFStates[i, :]  .= wf.States_WF[start_idx, :]
    end

    # Downwind step (same as basic)
    @inbounds @simd for i in 1:size(wf.States_OP, 1)
        buffers.step_dw[i] = sim.time_step * sim.dyn.advection * wf.States_WF[i, 1]
    end
    @inbounds @simd for i in 1:size(wf.States_OP, 1)
        wf.States_OP[i, 4] += buffers.step_dw[i]
    end

    # Crosswind deflection (reuse centerline!)
    centerline!(buffers.deflection, wf.States_OP, wf.States_T, wf.States_WF, floris, wf.D[1])
    @inbounds for i in 1:size(wf.States_OP, 1)
        buffers.step_cw[i, 1] = buffers.deflection[i, 1] - wf.States_OP[i, 5]
        buffers.step_cw[i, 2] = buffers.deflection[i, 2] - wf.States_OP[i, 6]
        wf.States_OP[i, 5] = buffers.deflection[i, 1]
        wf.States_OP[i, 6] = buffers.deflection[i, 2]
    end

    # Transform to world coordinates
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

    # Circular shift (manual) for OPs and restore turbine starts
    _circshift_and_restore!(wf.States_OP, buffers.tmpOPStates, wf.StartI, buffers.temp_states_op)
    _circshift_and_restore!(wf.States_T,  buffers.tmpTStates,  wf.StartI, buffers.temp_states_t)

    # Wind field circular shift with averaging weights
    wNew, wOld = sim.dyn.op_iter_weights  # Expect length 2 vector
    # Perform shift into temp buffer
    _circshift_and_restore!(wf.States_WF, buffers.tmpWFStates, wf.StartI, buffers.temp_states_wf)
    # Average: existing (already shifted & restored) is treated as "new" part
    @inbounds @simd for i in 1:size(wf.States_WF, 1)
        @inbounds for j in 1:size(wf.States_WF, 2)
            wf.States_WF[i, j] = wOld * buffers.temp_states_wf[i, j] + wNew * wf.States_WF[i, j]
        end
    end
    # Restore the turbine starting rows again (MATLAB restored after averaging)
    @inbounds for i in 1:size(wf.StartI, 2)
        start_idx = wf.StartI[1, i]
        @views wf.States_WF[start_idx, :] .= buffers.tmpWFStates[i, :]
    end

    # Reorder if downstream positions got unsorted
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

Check if observation points are in order and reorder them if not, using pre-allocated buffers.
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

# function T = iterateOPs(T,Sim,paramFLORIS,paramFLORIDyn)
# %ITERATEOPS propagates the OPs downstream and deletes the oldest ones.


# %% Save turbine OPs
# tmpOPStates = T.States_OP(T.StartI,:);
# tmpTStates  = T.States_T(T.StartI,:);
# tmpWFSTates = T.States_WF(T.StartI,:);
# %% Shift states
# % Calculate downwind step and apply to wake coordinate system
# step_dw = Sim.TimeStep * T.States_WF(:,1) * Sim.Dyn.Advection;
# T.States_OP(:,4) = T.States_OP(:,4) + step_dw;

# % Calculate crosswind step and apply to wake coordinate system
# %   TODO If turbines with different diameters are used, this has to be run
# %   individually (in parallel) for all turbines.
# deflection  = Centerline(...
#     T.States_OP,T.States_T,T.States_WF,paramFLORIS,T.D(1));
# step_cw     = (deflection - T.States_OP(:,5:6))*1;

# T.States_OP(:,5:6) = deflection*1;

# % Apply dw and cw step to the world coordinate system
# phiW = angSOWFA2world(T.States_WF(:,2));

# T.States_OP(:,1) = T.States_OP(:,1) + ...
#     cos(phiW) .* step_dw - ...
#     sin(phiW) .* step_cw(:,1);
# T.States_OP(:,2) = T.States_OP(:,2) + ...
#     sin(phiW) .* step_dw + ...
#     cos(phiW) .* step_cw(:,1);
# T.States_OP(:,3) = T.States_OP(:,3) + step_cw(:,2);

# %% Circshift & Init first OPs
# %   OPs
# T.States_OP = circshift(T.States_OP,1);
# T.States_OP(T.StartI,:) = tmpOPStates;
# %   Turbines
# T.States_T  = circshift(T.States_T,1);
# T.States_T(T.StartI,:)  = tmpTStates;
# %   Wind Field
# % Getting weights for averaging
# wNew = Sim.Dyn.OPiterWeights(1);
# wOld = Sim.Dyn.OPiterWeights(2);
# T.States_WF = wOld*T.States_WF + wNew*circshift(T.States_WF,1);
# T.States_WF(T.StartI,:) = tmpWFSTates;

# %% Check if OPs are in order
# for iT=1:T.nT
#     [~,indOP] = sort(T.States_OP(T.StartI(iT)+(0:(T.nOP-1)),4));
#     if ~issorted(indOP)
#         warning(['OPs overtaking, consider increasing the weight on ' ...
#             'the old wind field state in setup > OP / wind field propagation.'])
#         T.States_OP(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
#             T.States_OP(T.StartI(iT) + indOP - 1,:);
#         T.States_T(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
#             T.States_T(T.StartI(iT) + indOP - 1,:);
#         T.States_WF(T.StartI(iT)+(0:(T.nOP-1)),:) = ...
#             T.States_WF(T.StartI(iT) + indOP - 1,:);
#     end
# end

# end
