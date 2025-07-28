# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getMeasurements(mx::Matrix, my::Matrix, nM::I        # Update StartI to include the grid point
        GP.StartI = hcat(wf.StartI, [wf.StartI[end] + 1])
        # Update arrays to accommodate the additional grid point
        GP.States_OP = wf.States_OP
        GP.posBase = vcat(wf.posBase, reshape([xGP, yGP, 0.0], 1, 3))  # Add grid point position
        GP.nOP = wf.nOP
        GP.intOPs = interpolateOPs(GP)
        
        GP.posNac = vcat(wf.posNac, reshape([0.0, 0.0, zh], 1, 3))  # Add grid point nacelle position
        GP.States_WF = wf.States_WF
        GP.States_T = vcat(wf.States_T, zeros(1, size(wf.States_T, 2)))  # Add row for grid point
        GP.D = vcat(wf.D, [0.0])  # Add 0 diameter for grid point (it's not a real turbine), wf::WindFarm, set::Settings, floris::Floris, wind::Wind)

Wrapper function to disguise the grid points as turbines with one rotor
point and experience almost the same calculations as the rotor points in
the simulation.

Single thread version

# Arguments
- `mx::Matrix`: X-coordinates of grid points
- `my::Matrix`: Y-coordinates of grid points  
- `nM::Int`: Number of measurements
- `zh::Real`: Hub height
- `wf::WindFarm`: Wind farm object containing turbine data
- `set::Settings`: Settings object containing simulation parameters
- `floris::Floris`: FLORIS model parameters
- `wind::Wind`: Wind field configuration

# Returns
- `mz::Array{Float64,3}`: 3D array of measurements with dimensions (size(mx,1), size(mx,2), nM)
"""
function getMeasurements(mx, my, nM, zh, wf::WindFarm, set::Settings, floris::Floris, wind::Wind)
    size_mx = size(mx)
    mz = zeros(size_mx[1], size_mx[2], nM)
    # Single-threaded loop (can be parallelized with @threads or Distributed.@distributed)
    for iGP in 1:length(mx)
        xGP = mx[iGP]
        yGP = my[iGP]
        
        GPdep = collect(1:wf.nT)  # Convert range to vector

        # Create T equivalent to get interpolated OPs
        # Initialize GP as a new WindFarm with a single grid point
        GP = deepcopy(wf)  # Start with a copy to get the right structure
        
        # The grid point depends on all original turbines
        GPdep = collect(1:wf.nT)  
        
        # Set up GP as a wind farm with wf.nT+1 turbines (original turbines + 1 grid point)
        # The grid point is the last "turbine" and depends on all the original ones
        GP.nT = wf.nT + 1  # Original turbines + 1 grid point
        GP.dep = Vector{Vector{Int64}}(undef, GP.nT)
        
        # Initialize dependencies: original turbines are independent, grid point depends on all
        for i in 1:wf.nT
            GP.dep[i] = Int64[]  # Original turbines are independent (no dependencies)
        end
        
        # Grid point (last "turbine") depends on all original turbines
        GP.dep[end] = GPdep
        
        # Update StartI to include the grid point
        GP.StartI = hcat(wf.StartI, [wf.StartI[end] + 1])
        # Update arrays to accommodate the additional grid point
        GP.States_OP = wf.States_OP
        GP.posBase = vcat(wf.posBase, reshape([xGP, yGP, 0.0], 1, 3))  # Add grid point position
        GP.nOP = wf.nOP
        GP.intOPs = interpolateOPs(GP)
        
        GP.posNac = vcat(wf.posNac, reshape([0.0, 0.0, zh], 1, 3))  # Add grid point nacelle position
        GP.States_WF = wf.States_WF
        GP.States_T = vcat(wf.States_T, zeros(1, size(wf.States_T, 2)))  # Add row for grid point
        GP.D = vcat(wf.D, [0.0])  # Add 0 diameter for grid point (it's not a real turbine)

        tmpM, _ = setUpTmpWFAndRun(set, GP, floris, wind)
        
        # Extract only the result for the grid point (last "turbine")
        gridPointResult = tmpM[end, :]
        
        # Convert linear index to subscripts
        rw, cl = divrem(iGP - 1, size_mx[1])
        rw += 1
        cl += 1
        mz[rw, cl, 1:3] = gridPointResult
    end
    
    return mz
end

"""
    plotFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris)

Generate full flow field plot data by calculating measurements across a grid.

# Arguments
- `wf::WindFarm`: Wind farm object containing turbine data
- `wind::Wind`: Wind field configuration  
- `set::Settings`: Settings object containing simulation parameters
- `floris::Floris`: FLORIS model parameters

# Returns
- `Z::Array{Float64,3}`: 3D array of flow field measurements
- `X::Matrix{Float64}`: X-coordinate grid
- `Y::Matrix{Float64}`: Y-coordinate grid
"""
function plotFlowField(set::Settings, wf::WindFarm, wind::Wind, floris::Floris)
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
    println(size(X), size(Y))
    println("nM: ", nM, " zh: ", zh)
    Z = getMeasurements(X, Y, nM, zh, wf, set, floris, wind)
    
    return Z, X, Y
end
