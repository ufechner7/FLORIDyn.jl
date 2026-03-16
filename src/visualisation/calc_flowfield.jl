# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Base.Threads

# Forward declarations avoid early-binding world-age diagnostics in static analysis.
function create_thread_buffers end
function update_thread_buffers! end
function getMeasurements end
function calcFlowField end

"""
    ThreadBuffers

Thread-local buffers for parallel flow field computation.

# Fields
- `thread_buffers::Vector{WindFarm}`: Thread-local WindFarm objects for each thread
- `thread_unified_buffers::Vector{UnifiedBuffers}`: Thread-local UnifiedBuffers for computation
"""
struct ThreadBuffers
    thread_buffers::Vector{WindFarm}
    thread_unified_buffers::Vector{UnifiedBuffers}
end

# UnifiedBuffers struct and create_unified_buffers function are defined in floridyn_cl/structs.jl

function create_thread_buffers(wf::WindFarm, nth::Int, floris::Floris)
    # Pre-allocate unified buffers for each thread (contains GP)
    thread_unified_buffers = Vector{UnifiedBuffers}(undef, nth)
    thread_buffers = Vector{WindFarm}(undef, nth)

    for tid in 1:nth
        # Create unified buffers with proper FLORIS parameters for this thread
        ub = create_unified_buffers(wf, floris)
        thread_unified_buffers[tid] = ub
        # Use the prebuilt GP from unified buffers as the thread's WindFarm buffer
        thread_buffers[tid] = ub.gp
    end

    return ThreadBuffers(thread_buffers, thread_unified_buffers)
end

function create_thread_buffers(wf::WindFarm, nth::Int)
    # Pre-allocate unified buffers for each thread (contains GP)
    thread_unified_buffers = Vector{UnifiedBuffers}(undef, nth)
    thread_buffers = Vector{WindFarm}(undef, nth)

    for tid in 1:nth
        # Create unified buffers with default rotor discretization
        ub = create_unified_buffers(wf, 50)
        thread_unified_buffers[tid] = ub
        # Use the prebuilt GP from unified buffers as the thread's WindFarm buffer
        thread_buffers[tid] = ub.gp
    end

    return ThreadBuffers(thread_buffers, thread_unified_buffers)
end

@doc """
    create_thread_buffers(wf::WindFarm, nth::Int, floris::Floris) -> ThreadBuffers
    create_thread_buffers(wf::WindFarm, nth::Int) -> ThreadBuffers

Create thread-local buffers for parallel flow field computation.

These methods pre-allocate all necessary data structures for each thread to avoid
race conditions and memory allocations during the parallel computation loop.

# Arguments
- `wf::WindFarm`: Original wind farm object to use as template
- `nth::Int`: Number of threads to create buffers for
- `floris::Floris`: FLORIS parameters for creating proper FLORIS buffers

# Returns
- `ThreadBuffers`: Struct containing all thread-local buffers

# Performance Notes
- Each thread gets its own copy of the WindFarm structure
- Pre-allocates all arrays to minimize allocations during computation
- Sets up dependency structure for virtual turbines at grid points
""" create_thread_buffers

function update_thread_buffers!(buffers::ThreadBuffers, wf::WindFarm)
    # Update wind field states in all thread-local WindFarm objects
    for tid in 1:length(buffers.thread_buffers)
        thread_wf = buffers.thread_buffers[tid]
        
        # Update wind field states matrix (non-allocating in-place assignment)
        # Note: Only update the original turbine portion, not the extra grid point
        original_size = size(wf.States_WF)
        thread_wf.States_WF[1:original_size[1], 1:original_size[2]] .= wf.States_WF
        
        # Update optional interpolation coefficient matrices if they exist
        if hasfield(typeof(wf), :C_Vel) && isdefined(wf, :C_Vel)
            if hasfield(typeof(thread_wf), :C_Vel) && isdefined(thread_wf, :C_Vel)
                # Only update the original turbine portion
                original_nT = wf.nT
                thread_wf.C_Vel[1:original_nT, :] .= wf.C_Vel
            end
        end
        
        if hasfield(typeof(wf), :C_Dir) && isdefined(wf, :C_Dir)
            if hasfield(typeof(thread_wf), :C_Dir) && isdefined(thread_wf, :C_Dir)
                # Only update the original turbine portion
                original_nT = wf.nT
                thread_wf.C_Dir[1:original_nT, :] .= wf.C_Dir
            end
        end
    end
    
    return nothing
end

    @doc """
        update_thread_buffers!(buffers::ThreadBuffers, wf::WindFarm) -> Nothing

    Update wind field states in all thread-local wind farm buffers without allocating memory.

    This function efficiently updates the wind field states (`States_WF`) and optional interpolation
    coefficients (`C_Vel`, `C_Dir`) in all thread-local WindFarm objects to match the current
    wind conditions from the master WindFarm object.

    # Arguments
    - `buffers::ThreadBuffers`: Thread-local buffers containing WindFarm copies for each thread
    - `wf::WindFarm`: Master WindFarm object with updated wind field states

    # Performance Notes
    - Uses in-place assignment (`.=`) to avoid memory allocations
    - Only updates wind-related fields, preserving other thread-specific modifications
    - Updates all threads' buffers to maintain consistency
    - Handles optional interpolation coefficient matrices when present

    # See Also
    - `create_thread_buffers`: Create initial thread-local buffers
    - `getMeasurements`: Parallel flow field computation using thread buffers
    """ update_thread_buffers!



function getMeasurements(buffers, mx, my, nM, zh, wf::WindFarm, set::Settings, floris::Floris, wind::Wind)
    size_mx = size(mx)
    mz = zeros(size_mx[1], size_mx[2], nM)
    
    if length(buffers.thread_buffers) == 1
        # Single-threaded loop when only one buffer is provided
        for iGP in eachindex(mx)
            # Use the single available buffer
            GP = buffers.thread_buffers[1]
            unified_buffers = buffers.thread_unified_buffers[1]

            xGP = mx[iGP]
            yGP = my[iGP]

            GP.posBase[end, 1] = xGP
            GP.posBase[end, 2] = yGP
            GP.posBase[end, 3] = 0.0

            GP.posNac[end, 1] = 0.0
            GP.posNac[end, 2] = 0.0
            GP.posNac[end, 3] = zh

            # Reset the grid point state
            for j in axes(GP.States_T, 2)
                GP.States_T[end, j] = 0.0
            end

            interpolateOPs!(unified_buffers, GP.intOPs, GP)

            setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)
            tmpM = unified_buffers.M_buffer

            @views gridPointResult = tmpM[end, :]

            # Map linear index to (row, col) for column-major arrays
            q, r = divrem(iGP - 1, size_mx[1])
            row = r + 1
            col = q + 1
            @views mz[row, col, 1:3] .= gridPointResult
        end
    else
        # Parallel loop using @threads
        @threads :static for iGP in eachindex(mx)
            # Get thread-local buffers
            tid = threadid()
            GP = buffers.thread_buffers[tid]
            unified_buffers = buffers.thread_unified_buffers[tid]
            
            xGP = mx[iGP]
            yGP = my[iGP]
            
            # Thread-safe updates: create copies to avoid modifying shared arrays
            GP.posBase[end, 1] = xGP
            GP.posBase[end, 2] = yGP
            GP.posBase[end, 3] = 0.0
            
            GP.posNac[end, 1] = 0.0
            GP.posNac[end, 2] = 0.0
            GP.posNac[end, 3] = zh
            
            # Reset the grid point state (thread-safe element-wise assignment)
            for j in axes(GP.States_T, 2)
                GP.States_T[end, j] = 0.0
            end
            
            # Recalculate interpolated OPs for the updated geometry (non-allocating)
            interpolateOPs!(unified_buffers, GP.intOPs, GP)

            setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)
            tmpM = unified_buffers.M_buffer
            
            # Extract only the result for the grid point (last "turbine")
            @views gridPointResult = tmpM[end, :]
            
            # Convert linear index to subscripts (thread-safe)
            # Map linear index to (row, col) for column-major arrays
            q, r = divrem(iGP - 1, size_mx[1])
            row = r + 1
            col = q + 1
            
            # Thread-safe assignment using @views to avoid race conditions
            @views mz[row, col, 1:3] .= gridPointResult
        end
    end

    return mz
end

@doc """
    getMeasurements(buffers::ThreadBuffers, mx::Matrix, my::Matrix, nM::Int, zh::Real,
                    wf::WindFarm, set::Settings, floris::Floris, wind::Wind) -> Array{Float64,3}

Calculate flow field measurements at specified grid points by treating them as virtual turbines.

This function computes flow field properties (velocity reduction, added turbulence, effective wind speed)
at grid points by creating virtual turbines at each location and running the FLORIS wake model.
Each grid point is treated as a turbine that depends on all real turbines in the wind farm,
allowing wake effects to be captured in the flow field visualization.

# Arguments
- `buffers::ThreadBuffers`: Pre-allocated thread-local buffers created with [`create_thread_buffers`](@ref);
    for Julia 1.12 use `create_thread_buffers(wf, nthreads() + 1, floris)`; for single-thread use `create_thread_buffers(wf, 1, floris)`
- `mx::Matrix`: X-coordinates of grid points (m)
- `my::Matrix`: Y-coordinates of grid points (m)
- `nM::Int`: Number of measurements to compute (typically 3)
- `zh::Real`: Hub height for measurements (m)
- `wf::WindFarm`: Wind farm object containing turbine data. See: [`WindFarm`](@ref)
- `set::Settings`: Settings object containing simulation parameters. See: [`Settings`](@ref)
- `floris::Floris`: FLORIS model parameters for wake calculations. See: [`Floris`](@ref)
- `wind::Wind`: Wind field configuration. See: [`Wind`](@ref)

# Returns
- `mz::Array{Float64,3}`: 3D array of measurements with dimensions `(size(mx,1), size(mx,2), nM)`
  - `mz[:,:,1]`: Velocity reduction
  - `mz[:,:,2]`: Added turbulence intensity
  - `mz[:,:,3]`: Effective wind speed

# Performance Notes
- Multi-threaded implementation using `@threads` for parallel processing of grid points when more than one buffer is provided
- With a single buffer (`length(buffers.thread_buffers) == 1`), runs in a single-threaded loop
- Each grid point requires a full wind farm simulation, so computation time scales with grid size
- Uses thread-local buffers created by [`create_thread_buffers`](@ref) to avoid race conditions
- On Julia 1.12 create `nthreads() + 1` buffers to accommodate thread indexing

# See Also
- `calcFlowField`: Higher-level function that uses this to create complete flow field data
- [`setUpTmpWFAndRun!`](@ref): Underlying simulation function used for each grid point
""" getMeasurements

function calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris;
                       plt=nothing, vis=nothing)
    # Preallocate field
    nM = 3

    # Use vis struct fields if provided, otherwise fall back to defaults
    if vis !== nothing
        fieldLims = [vis.field_limits_min';
                     vis.field_limits_max']  # [xmin ymin zmin; xmax ymax zmax]
        fieldRes = vis.field_resolution
    else
        # Default values for backward compatibility
        fieldLims = [0.0 0.0 0.0;
                     3000.0 3000.0 400.0]  # [xmin ymin zmin; xmax ymax zmax]
        fieldRes = 20.0  # Resolution of the field in m
    end
    
    xAx = fieldLims[1,1]:fieldRes:fieldLims[2,1]
    yAx = fieldLims[1,2]:fieldRes:fieldLims[2,2]

    # Create coordinate grids (Julia equivalent of meshgrid)
    X = repeat(collect(xAx)', length(yAx), 1)
    Y = repeat(collect(yAx), 1, length(xAx))

    # Get hub height from first turbine
    zh = wf.posNac[1, 3]
    
    # Get data
    if set.threading
        # Disable garbage collection before multithreading if not parallel
        if ! set.parallel
            GC.enable(false)
        end
        try
            # Create thread-local buffers using the new function with FLORIS parameters
            buffers = create_thread_buffers(wf, nthreads() + 1, floris)
            Z = getMeasurements(buffers, X, Y, nM, zh, wf, set, floris, wind)
        finally
            # Re-enable garbage collection after multithreading if not parallel
            if ! set.parallel
                GC.enable(true)
            end
        end
    else
        # Single-threaded path using getMeasurements with a single buffer
        buffers = create_thread_buffers(wf, 1, floris)
        Z = getMeasurements(buffers, X, Y, nM, zh, wf, set, floris, wind)
    end

    return Z, X, Y
end

@doc """
    calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris; plt=nothing)

Generate full flow field plot data by calculating measurements across a grid.

This function creates a rectangular grid over the wind farm domain and calculates flow field
properties at each grid point by treating them as virtual turbines. The computation can be
performed in parallel if `set.threading` is true.

# Arguments
- `set::Settings`: Settings object containing simulation parameters
- `wf::WindFarm`: Wind farm object containing turbine data
- `wind::Wind`: Wind field configuration
- `floris::Floris`: FLORIS model parameters

# Keyword Arguments
- `plt=nothing`: Plot object for garbage collection control.
- `vis=nothing`: Visualization configuration object containing field limits and resolution settings.

# Returns
- `Z::Array{Float64,3}`: 3D array of flow field measurements with dimensions `(ny, nx, 3)`
- `X::Matrix{Float64}`: X-coordinate grid (m)
- `Y::Matrix{Float64}`: Y-coordinate grid (m)

# See Also
- `getMeasurements`: Function used internally to compute the flow field
- [`plotFlowField`](@ref): Visualization function for the generated data
""" calcFlowField
