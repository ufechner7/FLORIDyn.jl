# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Base.Threads

"""
    ThreadBuffers

Thread-local buffers for parallel flow field computation.

# Fields
- `thread_buffers::Vector{WindFarm}`: Thread-local WindFarm objects for each thread
- `thread_comp_buffers::Vector{NamedTuple}`: Computation buffers containing distance and sorting arrays
"""
struct ThreadBuffers
    thread_buffers::Vector{WindFarm}
    thread_comp_buffers::Vector{@NamedTuple{
        dist_buffer::Vector{Float64}, 
        sorted_indices_buffer::Vector{Int},
        M_buffer::Matrix{Float64},
        iTWFState_buffer::Vector{Float64},
        tmp_Tpos_buffer::Matrix{Float64},
        tmp_WF_buffer::Matrix{Float64},
        tmp_Tst_buffer::Matrix{Float64},
        dists_buffer::Vector{Float64},
        plot_WF_buffer::Matrix{Float64},
        plot_OP_buffer::Matrix{Float64}
    }}
end

# UnifiedBuffers struct and create_unified_buffers function are defined in floridyn_cl/structs.jl

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
    # Create thread-local buffers to avoid race conditions
    thread_buffers = [deepcopy(wf) for _ in 1:nth]
    
    # Pre-allocate arrays for each thread buffer
    original_nT = wf.nT
    
    # Pre-allocate computation buffers for each thread using unified buffers
    thread_comp_buffers = Vector{@NamedTuple{
        dist_buffer::Vector{Float64}, 
        sorted_indices_buffer::Vector{Int},
        M_buffer::Matrix{Float64},
        iTWFState_buffer::Vector{Float64},
        tmp_Tpos_buffer::Matrix{Float64},
        tmp_WF_buffer::Matrix{Float64},
        tmp_Tst_buffer::Matrix{Float64},
        dists_buffer::Vector{Float64},
        plot_WF_buffer::Matrix{Float64},
        plot_OP_buffer::Matrix{Float64}
    }}(undef, nth)
    
    for tid in 1:nth
        GP = thread_buffers[tid]
        GP.nT = original_nT + 1  # Original turbines + 1 grid point
        
        # Pre-allocate dependency structure
        GP.dep = Vector{Vector{Int64}}(undef, GP.nT)
        for i in 1:original_nT
            GP.dep[i] = Int64[]  # Original turbines are independent (no dependencies)
        end
        GP.dep[end] = collect(1:original_nT)  # Grid point depends on all original turbines
        
        # Pre-allocate extended arrays
        GP.StartI = hcat(wf.StartI, [wf.StartI[end] + 1])
        GP.posBase = vcat(wf.posBase, zeros(1, 3))  # Will be updated for each grid point
        GP.posNac = vcat(wf.posNac, zeros(1, 3))   # Will be updated for each grid point
        GP.States_T = vcat(wf.States_T, zeros(1, size(wf.States_T, 2)))
        GP.D = vcat(wf.D, [0.0])  # Grid point has 0 diameter
        
        # Ensure all necessary fields are copied to avoid shared references
        if hasfield(typeof(wf), :intOPs)
            GP.intOPs = deepcopy(wf.intOPs)
        end
        
        # Pre-allocate buffers for non-allocating interpolateOPs!
        GP.intOPs = [zeros(length(GP.dep[iT]), 4) for iT in 1:GP.nT]
        
        # Create unified buffers for this thread
        unified_buffers = create_unified_buffers(wf)
        
        thread_comp_buffers[tid] = (
            dist_buffer = unified_buffers.dist_buffer,
            sorted_indices_buffer = unified_buffers.sorted_indices_buffer,
            M_buffer = unified_buffers.M_buffer,
            iTWFState_buffer = unified_buffers.iTWFState_buffer,
            tmp_Tpos_buffer = unified_buffers.tmp_Tpos_buffer,
            tmp_WF_buffer = unified_buffers.tmp_WF_buffer,
            tmp_Tst_buffer = unified_buffers.tmp_Tst_buffer,
            dists_buffer = unified_buffers.dists_buffer,
            plot_WF_buffer = unified_buffers.plot_WF_buffer,
            plot_OP_buffer = unified_buffers.plot_OP_buffer
        )
    end
    
    return ThreadBuffers(thread_buffers, thread_comp_buffers)
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
- [`setUpTmpWFAndRun`](@ref): Underlying simulation function used for each grid point
"""
function getMeasurements(mx, my, nM, zh, wf::WindFarm, set::Settings, floris::Floris, wind::Wind)
    size_mx = size(mx)
    mz = zeros(size_mx[1], size_mx[2], nM)
    
    # Pre-allocate GP buffer once outside the loop
    GP = deepcopy(wf)  # Create the buffer once
    
    # Pre-allocate arrays that will be modified for each grid point
    # These will be reused and updated for each iteration
    original_nT = wf.nT
    GP.nT = original_nT + 1  # Original turbines + 1 grid point
    
    # Pre-allocate dependency structure
    GP.dep = Vector{Vector{Int64}}(undef, GP.nT)
    for i in 1:original_nT
        GP.dep[i] = Int64[]  # Original turbines are independent (no dependencies)
    end
    GP.dep[end] = collect(1:original_nT)  # Grid point depends on all original turbines
    
    # Pre-allocate extended arrays
    GP.StartI = hcat(wf.StartI, [wf.StartI[end] + 1])
    GP.posBase = vcat(wf.posBase, zeros(1, 3))  # Will be updated for each grid point
    GP.posNac = vcat(wf.posNac, zeros(1, 3))   # Will be updated for each grid point
    GP.States_T = vcat(wf.States_T, zeros(1, size(wf.States_T, 2)))
    GP.D = vcat(wf.D, [0.0])  # Grid point has 0 diameter
    
    # Pre-allocate buffers for non-allocating interpolateOPs!
    GP.intOPs = [zeros(length(GP.dep[iT]), 4) for iT in 1:GP.nT]
    
    # Create a single unified buffer struct containing all arrays
    unified_buffers = create_unified_buffers(wf)

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
        interpolateOPs!(GP.intOPs, GP, unified_buffers)

        tmpM, _ = setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)
        
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
- [`setUpTmpWFAndRun`](@ref): Underlying simulation function used for each grid point
"""
function getMeasurementsP(mx, my, nM, zh, wf::WindFarm, set::Settings, floris::Floris, wind::Wind; alloc=nothing)
    size_mx = size(mx)
    a = @allocated mz = zeros(size_mx[1], size_mx[2], nM)
    if ! isnothing(alloc)
        alloc.gmp_mx += a
    end

    # for Julia 1.12 we need one extra thread
    nth = nthreads() + 1
    
    # Create thread-local buffers using the new function
    a = @allocated buffers = create_thread_buffers(wf, nth)
    if ! isnothing(alloc)
        alloc.gmp_buffers += a
    end
    gmp_alloc2 = Threads.Atomic{Int64}(0)
    
    # Parallel loop using @threads
    @threads :static for iGP in 1:length(mx)
        # Get thread-local buffers
        tid = threadid()
        GP = buffers.thread_buffers[tid]
        comp_buffers = buffers.thread_comp_buffers[tid]
        
        # Create a single unified buffer struct from the named tuple
        unified_buffers = UnifiedBuffers(
            comp_buffers.dist_buffer,
            comp_buffers.sorted_indices_buffer,
            comp_buffers.M_buffer,
            comp_buffers.iTWFState_buffer,
            comp_buffers.tmp_Tpos_buffer,
            comp_buffers.tmp_WF_buffer,
            comp_buffers.tmp_Tst_buffer,
            comp_buffers.dists_buffer,
            comp_buffers.plot_WF_buffer,
            comp_buffers.plot_OP_buffer
        )
        
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
        interpolateOPs!(GP.intOPs, GP, unified_buffers)

        a = @allocated tmpM, _ = setUpTmpWFAndRun!(unified_buffers, GP, set, floris, wind)

        Threads.atomic_add!(gmp_alloc2, a)
        
        # Extract only the result for the grid point (last "turbine")
        @views gridPointResult = tmpM[end, :]
        
        # Convert linear index to subscripts (thread-safe)
        rw, cl = divrem(iGP - 1, size_mx[1])
        rw += 1
        cl += 1
        
        # Thread-safe assignment using @views to avoid race conditions
        @views mz[rw, cl, 1:3] .= gridPointResult
    end
    if ! isnothing(alloc)
        alloc.gmp_alloc2 += gmp_alloc2[]  # Add total allocation to the main allocator
    end

    return mz
end

"""
    calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris; plt=nothing)

Generate full flow field plot data by calculating measurements across a grid.

This function creates a rectangular grid over the wind farm domain and calculates flow field
properties at each grid point by treating them as virtual turbines. The computation can be
performed in parallel if `set.parallel` is true, with automatic garbage collection management
for thread safety.

# Arguments
- `set::Settings`: Settings object containing simulation parameters
  - `set.parallel`: If true, uses parallel computation with multiple threads
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
# Calculate flow field with parallel computation and GC control
set.parallel = true
Z, X, Y = calcFlowField(set, wf, wind, floris; plt)

# Extract velocity reduction field
velocity_reduction = Z[:, :, 1]

# Extract effective wind speed field  
wind_speed = Z[:, :, 3]
```

# See Also
- [`getMeasurements`](@ref): Single-threaded implementation
- [`getMeasurementsP`](@ref): Multi-threaded implementation
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
    a = @allocated X = repeat(collect(xAx)', length(yAx), 1)
    b = @allocated Y = repeat(collect(yAx), 1, length(xAx))
    if ! isnothing(alloc)
        alloc.cff_X += a
        alloc.cff_Y += b
    end
    
    # Get hub height from first turbine
    zh = wf.posNac[1, 3]
    
    # Get data
    if set.threading
        # Disable garbage collection before multithreading if not parallel
        if ! set.parallel
            GC.enable(false)
        end
        try
            a = @allocated Z = getMeasurementsP(X, Y, nM, zh, wf, set, floris, wind; alloc)
            if ! isnothing(alloc)
                alloc.getMeasurementsP += a
            end
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
