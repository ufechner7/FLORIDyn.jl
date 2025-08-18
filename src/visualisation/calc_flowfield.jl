# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Base.Threads

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

"""
    create_thread_buffers(wf::WindFarm, nth::Int, floris::Floris) -> ThreadBuffers

Create thread-local buffers for parallel flow field computation with FLORIS parameters.

This function pre-allocates all necessary data structures for each thread to avoid
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
"""
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

"""
    create_thread_buffers(wf::WindFarm, nth::Int) -> ThreadBuffers

Create thread-local buffers for parallel flow field computation.

This function pre-allocates all necessary data structures for each thread to avoid
race conditions and memory allocations during the parallel computation loop.

# Arguments
- `wf::WindFarm`: Original wind farm object to use as template
- `nth::Int`: Number of threads to create buffers for

# Returns
- `ThreadBuffers`: Struct containing all thread-local buffers

# Performance Notes
- Each thread gets its own copy of the WindFarm structure
- Pre-allocates all arrays to minimize allocations during computation
- Sets up dependency structure for virtual turbines at grid points
"""
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

"""
    update_thread_buffers!(buffers::ThreadBuffers, wf::WindFarm) -> Nothing

Update wind field states in all thread-local wind farm buffers without allocating memory.

This function efficiently updates the wind field states (`States_WF`) and optional interpolation 
coefficients (`C_Vel`, `C_Dir`) in all thread-local WindFarm objects to match the current 
wind conditions from the master WindFarm object. This is useful when wind conditions change 
during simulation and the thread buffers need to be synchronized.

# Arguments
- `buffers::ThreadBuffers`: Thread-local buffers containing WindFarm copies for each thread
- `wf::WindFarm`: Master WindFarm object with updated wind field states

# Performance Notes
- Uses in-place assignment (`.=`) to avoid memory allocations
- Only updates wind-related fields, preserving other thread-specific modifications
- Updates all threads' buffers to maintain consistency
- Handles optional interpolation coefficient matrices when present

# Fields Updated
- `States_WF`: Wind field states matrix (velocity, direction, turbulence intensity)
- `C_Vel`: Velocity interpolation coefficients (if present in WindFarm type)  
- `C_Dir`: Direction interpolation coefficients (if present in WindFarm type)

# Example
```julia
# Create thread buffers
buffers = create_thread_buffers(wf, nthreads())

# ... wind conditions change ...

# Update all thread buffers with new wind states (non-allocating)
update_thread_buffers!(buffers, wf)

# Continue with flow field computation using updated buffers
Z = getMeasurementsP(buffers, X, Y, nM, zh, wf, set, floris, wind)
```

# See Also
- [`create_thread_buffers`](@ref): Create initial thread-local buffers
- [`getMeasurementsP`](@ref): Parallel flow field computation using thread buffers
"""
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

"""
    getMeasurements(mx::Matrix, my::Matrix, nM::Int, zh::Real, wf::WindFarm, set::Settings, 
                    floris::Floris, wind::Wind) -> Array{Float64,3}

Calculate flow field measurements at specified grid points by treating them as virtual turbines.

This function computes flow field properties (velocity reduction, added turbulence, effective wind speed)
at grid points by creating virtual turbines at each location and running the FLORIS wake model.
Each grid point is treated as a turbine that depends on all real turbines in the wind farm,
allowing wake effects to be captured in the flow field visualization.

# Arguments
- `mx::Matrix`: X-coordinates of grid points (m)
- `my::Matrix`: Y-coordinates of grid points (m)  
- `nM::Int`: Number of measurements to compute (typically 3)
- `zh::Real`: Hub height for measurements (m)
- `wf::WindFarm`: Wind farm object containing turbine data
  - `wf.nT`: Number of real turbines
  - `wf.StartI`: Starting indices for turbine data
  - `wf.posBase`, `wf.posNac`: Turbine positions
  - `wf.States_*`: Turbine state matrices
- `set::Settings`: Settings object containing simulation parameters
- `floris::Floris`: FLORIS model parameters for wake calculations
- `wind::Wind`: Wind field configuration

# Returns
- `mz::Array{Float64,3}`: 3D array of measurements with dimensions `(size(mx,1), size(mx,2), nM)`
  - `mz[:,:,1]`: Velocity reduction
  - `mz[:,:,2]`: Added turbulence intensity
  - `mz[:,:,3]`: Effective wind speed

# Algorithm
For each grid point:
1. Creates a temporary wind farm with all original turbines plus one virtual turbine at the grid point
2. Sets the virtual turbine to depend on all real turbines (to capture wake effects)
3. Runs the FLORIS simulation to compute wake-affected flow properties
4. Extracts the result for the virtual turbine position

# Performance Notes
- Single-threaded implementation (can be parallelized with `@threads` for large grids)
- Each grid point requires a full wind farm simulation, so computation time scales with grid size
- Uses pre-allocated buffer to avoid repeated `deepcopy` operations for better performance

# Example
```julia
# Create a 10x10 grid from 0 to 1000m
x_range = 0:100:1000
y_range = 0:100:1000
mx = repeat(collect(x_range)', length(y_range), 1)
my = repeat(collect(y_range), 1, length(x_range))

# Calculate 3 measurements at 90m hub height
mz = getMeasurements(mx, my, 3, 90.0, wind_farm, settings, floris_model, wind_config)

# Extract effective wind speed field
wind_speed_field = mz[:, :, 3]
```

# See Also
- [`calcFlowField`](@ref): Higher-level function that uses this to create complete flow field data
- [`setUpTmpWFAndRun!`](@ref): Underlying simulation function used for each grid point
"""
function getMeasurements(mx, my, nM, zh, wf::WindFarm, set::Settings, floris::Floris, wind::Wind; buffers=nothing)
    size_mx = size(mx)
    mz = zeros(size_mx[1], size_mx[2], nM)
    
    # Create a single unified buffer struct containing all arrays, including GP buffer
    unified_buffers = create_unified_buffers(wf)
    GP = unified_buffers.gp

    # Single-threaded loop (can be parallelized with @threads or DistributedNext.@distributed)
    for iGP in 1:length(mx)
        xGP = mx[iGP]
        yGP = my[iGP]
        
        # Update only the grid point position in the pre-allocated arrays
        GP.posBase[end, :] = [xGP, yGP, 0.0]  # Update grid point base position
        GP.posNac[end, :] = [0.0, 0.0, zh]    # Update grid point nacelle position
        
        # Reset the grid point state (other turbine states remain unchanged)
        GP.States_T[end, :] .= 0.0
        
        # Recalculate interpolated OPs for the updated geometry (non-allocating)
    interpolateOPs!(unified_buffers, GP.intOPs, GP)

        setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)
        tmpM = unified_buffers.M_buffer
        
        # Extract only the result for the grid point (last "turbine")
        @views gridPointResult = tmpM[end, :]
        
        # Convert linear index to subscripts
        rw, cl = divrem(iGP - 1, size_mx[1])
        rw += 1
        cl += 1
        mz[rw, cl, 1:3] = gridPointResult
    end
    
    return mz
end


"""
    getMeasurementsP(mx::Matrix, my::Matrix, nM::Int, zh::Real, wf::WindFarm, set::Settings, 
                     floris::Floris, wind::Wind) -> Array{Float64,3}

Calculate flow field measurements at specified grid points by treating them as virtual turbines.

This function computes flow field properties (velocity reduction, added turbulence, effective wind speed)
at grid points by creating virtual turbines at each location and running the FLORIS wake model.
Each grid point is treated as a turbine that depends on all real turbines in the wind farm,
allowing wake effects to be captured in the flow field visualization.

# Arguments
- `mx::Matrix`: X-coordinates of grid points (m)
- `my::Matrix`: Y-coordinates of grid points (m)  
- `nM::Int`: Number of measurements to compute (typically 3)
- `zh::Real`: Hub height for measurements (m)
- `wf::WindFarm`: Wind farm object containing turbine data
  - `wf.nT`: Number of real turbines
  - `wf.StartI`: Starting indices for turbine data
  - `wf.posBase`, `wf.posNac`: Turbine positions
  - `wf.States_*`: Turbine state matrices
- `set::Settings`: Settings object containing simulation parameters
- `floris::Floris`: FLORIS model parameters for wake calculations
- `wind::Wind`: Wind field configuration

# Returns
- `mz::Array{Float64,3}`: 3D array of measurements with dimensions `(size(mx,1), size(mx,2), nM)`
  - `mz[:,:,1]`: Velocity reduction
  - `mz[:,:,2]`: Added turbulence intensity
  - `mz[:,:,3]`: Effective wind speed

# Algorithm
For each grid point:
1. Creates a temporary wind farm with all original turbines plus one virtual turbine at the grid point
2. Sets the virtual turbine to depend on all real turbines (to capture wake effects)
3. Runs the FLORIS simulation to compute wake-affected flow properties
4. Extracts the result for the virtual turbine position

# Performance Notes
- Multi-threaded implementation using `@threads` for parallel processing of grid points
- Each grid point requires a full wind farm simulation, so computation time scales with grid size
- Uses thread-local buffers created by [`create_thread_buffers`](@ref) to avoid race conditions
- Scales well with the number of available CPU cores

# Example
```julia
# Create a 10x10 grid from 0 to 1000m
x_range = 0:100:1000
y_range = 0:100:1000
mx = repeat(collect(x_range)', length(y_range), 1)
my = repeat(collect(y_range), 1, length(x_range))

# Calculate 3 measurements at 90m hub height
mz = getMeasurements(mx, my, 3, 90.0, wind_farm, settings, floris_model, wind_config)

# Extract effective wind speed field
wind_speed_field = mz[:, :, 3]
```

# See Also
- [`calcFlowField`](@ref): Higher-level function that uses this to create complete flow field data
- [`setUpTmpWFAndRun!`](@ref): Underlying simulation function used for each grid point
"""
function getMeasurementsP(buffers, mx, my, nM, zh, wf::WindFarm, set::Settings, floris::Floris, wind::Wind; alloc=nothing)
    size_mx = size(mx)
    mz = zeros(size_mx[1], size_mx[2], nM)
       
    # Parallel loop using @threads
    @threads :static for iGP in 1:length(mx)
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
        for j in 1:size(GP.States_T, 2)
            GP.States_T[end, j] = 0.0
        end
        
        # Recalculate interpolated OPs for the updated geometry (non-allocating)
        interpolateOPs!(unified_buffers, GP.intOPs, GP)

        setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)
        tmpM = unified_buffers.M_buffer
        
        # Extract only the result for the grid point (last "turbine")
        @views gridPointResult = tmpM[end, :]
        
        # Convert linear index to subscripts (thread-safe)
        rw, cl = divrem(iGP - 1, size_mx[1])
        rw += 1
        cl += 1
        
        # Thread-safe assignment using @views to avoid race conditions
        @views mz[rw, cl, 1:3] .= gridPointResult
    end

    return mz
end

"""
    calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris; plt=nothing, alloc=nothing)

Generate full flow field plot data by calculating measurements across a grid.

This function creates a rectangular grid over the wind farm domain and calculates flow field
properties at each grid point by treating them as virtual turbines. The computation can be
performed in parallel if `set.threading` is true, with automatic garbage collection management
for thread safety.

# Arguments
- `set::Settings`: Settings object containing simulation parameters
  - `set.threading`: If true, uses multi-threaded computation with `@threads`
  - `set.parallel`: If true, enables parallel-specific optimizations
- `wf::WindFarm`: Wind farm object containing turbine data
- `wind::Wind`: Wind field configuration  
- `floris::Floris`: FLORIS model parameters

# Keyword Arguments
- `plt=nothing`: Plot object for garbage collection control. If provided and `set.parallel` is true,
  automatically calls `plt.GC.enable(false)` before multithreading and `plt.GC.enable(true)` 
  after completion to prevent PyCall-related segmentation faults during parallel execution.

# Returns
- `Z::Array{Float64,3}`: 3D array of flow field measurements with dimensions `(ny, nx, 3)`
  - `Z[:,:,1]`: Velocity reduction factor
  - `Z[:,:,2]`: Added turbulence intensity  
  - `Z[:,:,3]`: Effective wind speed (m/s)
- `X::Matrix{Float64}`: X-coordinate grid (m)
- `Y::Matrix{Float64}`: Y-coordinate grid (m)

# Performance Notes
- Uses parallel computation when `set.parallel=true` for significant speedup on multi-core systems
- Automatic garbage collection management prevents threading-related crashes
- Grid resolution is fixed at 20m with domain from [0,0] to [3000,3000] meters
- Hub height is taken from the first turbine in the wind farm

# Example
```julia
# Calculate flow field with threading and GC control
set.threading = true
Z, X, Y = calcFlowField(set, wf, wind, floris; plt)

# Extract velocity reduction field
velocity_reduction = Z[:, :, 1]

# Extract effective wind speed field  
wind_speed = Z[:, :, 3]
```

# See Also
- [`getMeasurements`](@ref): Unified function for single and multi-threaded flow field computation
- [`plotFlowField`](@ref): Visualization function for the generated data
"""
function calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris; plt=nothing, alloc=nothing)
    # Preallocate field
    nM = 3
    fieldLims = [0.0 0.0 0.0;
                 3000.0 3000.0 400.0]  # [xmin ymin zmin; xmax ymax zmax]
    fieldRes = 20.0  # Resolution of the field in m
    
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
            Z = getMeasurementsP(buffers, X, Y, nM, zh, wf, set, floris, wind)
        finally
            # Re-enable garbage collection after multithreading if not parallel
            if ! set.parallel
                GC.enable(true)
            end
        end
    else
        Z = getMeasurements(X, Y, nM, zh, wf, set, floris, wind)
    end
    
    return Z, X, Y
end
