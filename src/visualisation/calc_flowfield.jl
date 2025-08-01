# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using OhMyThreads

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
    dist_buffer = zeros(wf.nOP)
    sorted_indices_buffer = zeros(Int, wf.nOP)

    # Single-threaded loop (can be parallelized with @threads or Distributed.@distributed)
    for iGP in 1:length(mx)
        xGP = mx[iGP]
        yGP = my[iGP]
        
        # Update only the grid point position in the pre-allocated arrays
        GP.posBase[end, :] = [xGP, yGP, 0.0]  # Update grid point base position
        GP.posNac[end, :] = [0.0, 0.0, zh]    # Update grid point nacelle position
        
        # Reset the grid point state (other turbine states remain unchanged)
        GP.States_T[end, :] .= 0.0
        
        # Recalculate interpolated OPs for the updated geometry (non-allocating)
        interpolateOPs!(GP.intOPs, GP, dist_buffer, sorted_indices_buffer)

        tmpM, _ = setUpTmpWFAndRun(set, GP, floris, wind)
        
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
- Multi-threaded implementation using `OhMyThreads.jl` for efficient parallel processing
- Uses dynamic task scheduling (`scheduler=:dynamic`) for optimal load balancing across CPU cores
- Task-local memory allocation eliminates race conditions and reduces memory contention
- Each grid point requires a full wind farm simulation, so computation time scales with grid size
- Scales well with the number of available CPU cores and provides better performance than `@threads`

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
function getMeasurementsP(mx, my, nM, zh, wf::WindFarm, set::Settings, floris::Floris, wind::Wind)
    size_mx = size(mx)
    mz = zeros(size_mx[1], size_mx[2], nM)
    
    # Use OhMyThreads.tmap with dynamic scheduling for optimal load balancing
    original_nT = wf.nT
    
    # Process all grid points in parallel using tmap
    results = tmap(1:length(mx); scheduler=:dynamic) do iGP
        # Create task-local wind farm copy (each task gets its own copy)
        GP = deepcopy(wf)
        
        # Setup grid point as virtual turbine
        GP.nT = original_nT + 1
        
        # Pre-allocate dependency structure (task-local)
        GP.dep = Vector{Vector{Int64}}(undef, GP.nT)
        for i in 1:original_nT
            GP.dep[i] = Int64[]  # Original turbines are independent
        end
        GP.dep[end] = collect(1:original_nT)  # Grid point depends on all original turbines
        
        # Pre-allocate extended arrays (task-local)
        GP.StartI = hcat(wf.StartI, [wf.StartI[end] + 1])
        GP.posBase = vcat(wf.posBase, zeros(1, 3))
        GP.posNac = vcat(wf.posNac, zeros(1, 3))
        GP.States_T = vcat(wf.States_T, zeros(1, size(wf.States_T, 2)))
        GP.D = vcat(wf.D, [0.0])  # Grid point has 0 diameter
        
        # Pre-allocate interpolation buffers (task-local)
        GP.intOPs = [zeros(length(GP.dep[iT]), 4) for iT in 1:GP.nT]
        dist_buffer = zeros(wf.nOP)
        sorted_indices_buffer = zeros(Int, wf.nOP)
        
        # Get grid point coordinates
        xGP = mx[iGP]
        yGP = my[iGP]
        
        # Update grid point position
        GP.posBase[end, :] = [xGP, yGP, 0.0]
        GP.posNac[end, :] = [0.0, 0.0, zh]
        
        # Reset grid point state
        GP.States_T[end, :] .= 0.0
        
        # Recalculate interpolated OPs for updated geometry
        interpolateOPs!(GP.intOPs, GP, dist_buffer, sorted_indices_buffer)
        
        # Run simulation
        tmpM, _ = setUpTmpWFAndRun(set, GP, floris, wind)
        
        # Extract grid point result and compute indices
        gridPointResult = tmpM[end, :]
        rw, cl = divrem(iGP - 1, size_mx[1])
        rw += 1
        cl += 1
        
        # Return result with indices for thread-safe assignment
        (rw, cl, gridPointResult)
    end
    
    # Collect results into output array (serial section)
    for (rw, cl, result) in results
        mz[rw, cl, 1:3] = result
    end
    
    return mz
end

"""
    calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris)

Generate full flow field plot data by calculating measurements across a grid.

# Arguments
- `set::Settings`: Settings object containing simulation parameters
- `wf::WindFarm`: Wind farm object containing turbine data
- `wind::Wind`: Wind field configuration  
- `floris::Floris`: FLORIS model parameters

# Returns
- `Z::Array{Float64,3}`: 3D array of flow field measurements
- `X::Matrix{Float64}`: X-coordinate grid
- `Y::Matrix{Float64}`: Y-coordinate grid
"""
function calcFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris)
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
    if set.parallel
        Z = getMeasurementsP(X, Y, nM, zh, wf, set, floris, wind)
    else
        Z = getMeasurements(X, Y, nM, zh, wf, set, floris, wind)
    end
    
    return Z, X, Y
end
