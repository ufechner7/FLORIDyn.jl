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
    Z = getMeasurements(X, Y, nM, zh, wf, set, floris, wind)
    
    return Z, X, Y
end

"""
    plotFlowField(mx::Matrix, my::Matrix, mz::Array{Float64,3}; msr=1, title="Flow Field")

Plot a 2D contour of the flow field data.

# Arguments
- `mx::Matrix`: X-coordinate grid
- `my::Matrix`: Y-coordinate grid  
- `mz::Array{Float64,3}`: 3D array of measurements with dimensions (rows, cols, nM)
- `msr::Int`: Which measurement to plot (1, 2, or 3). Default is 3.

# Returns
- `nothing`

# Note
The measurement indices correspond to:
- 1: Velocity reduction
- 2: Added turbulence  
- 3: Effective wind speed

This function requires a plotting package like PyPlot.jl to be loaded and available as `plt`.
"""
function plotFlowField(plt, wf, mx, my, mz; msr=3, unit_test=false)
    # Extract the 2D slice for the specified measurement
    if msr > size(mz, 3)
        error("msr ($msr) exceeds number of measurements ($(size(mz, 3)))")
    end
    
    # Get the 2D slice
    mz_2d = mz[:, :, msr]
    
    # Try to use PyPlot if available
    try
        # This will work if PyPlot is loaded and plt is available
        if msr == 1
            figure = "Velocity Reduction"
        elseif msr == 2
            figure = "Added Turbulence"
        elseif msr == 3
            figure = "Effective Wind Speed"
        end
        title = figure
        n = 0.84
        fig = plt.figure(figure, figsize=(7.25n, 6n))
        vmin = 2.0; vmax = 10.0; n=40
        levels = range(vmin, stop=vmax, length=n+1)
        contour_plot = plt.contourf(my, mx, mz_2d, n; levels, cmap="inferno") # 40 levels, no lines
        #plt.axis("equal")
        cb = plt.colorbar()
        cb[:set_label](L"Wind speed~[ms^{-1}]", labelpad=3)

        # Plot the turbine rotors as short, thick lines (as seen from above)
        for i_T in 1:length(wf.D)
            # Compute yaw angle
            yaw = angSOWFA2world(wf.States_WF[wf.StartI[i_T], 2] - wf.States_T[wf.StartI[i_T], 2])

            # Rotation matrix
            R = [cos(yaw) -sin(yaw);
                sin(yaw) cos(yaw)]

            # Define rotor line endpoints before rotation (z ignored here)
            # Two points: (0, D/2) and (0, -D/2) along the vertical axis in local coords
            rotor_points = [0 0;
                            wf.D[i_T]/2 -wf.D[i_T]/2]

            # Apply rotation
            rot_pos = R * rotor_points

            # Add base position coordinates (broadcast)
            base_pos = wf.posBase[i_T, 1:2]  # (x,y)
            rot_pos .+= base_pos

            # Plot in 3D at height z=20 for both points
            ax = plt.gca()
            ax.plot(rot_pos[1, :], rot_pos[2, :], [20, 20], color="k", linewidth=3)
        end
        # Plot the OPs
        # Plot all points with size 5 and white filled marker
        n = 20
        plt.scatter(wf.States_OP[:, 1], wf.States_OP[:, 2], s=2, color="white", marker="o")
        plt.xlim(minimum(mx), maximum(mx))
        plt.ylim(minimum(mx), maximum(mx))

        # Plot every 10th point with size 15 and white filled marker
        plt.scatter(wf.States_OP[1:10:end, 1], wf.States_OP[1:10:end, 2], s=6, color="white", marker="o")
        
        plt.xlim(minimum(mx), maximum(mx))
        plt.ylim(minimum(mx), maximum(mx))
        plt.title(title)
        plt.xlabel("West-East [m]")
        plt.ylabel("South-North [m]")
        plt.tight_layout()
        println("Contour plot created successfully")
        if unit_test
            plt.pause(2)
            plt.close(fig)
        end
    catch e
        if isa(e, UndefVarError) && e.var == :plt
            @warn "PyPlot not available. Please load PyPlot.jl first with: using PyPlot; const plt = PyPlot"
            println("Data shape: ", size(mz_2d))
            println("X range: [", minimum(mx), ", ", maximum(mx), "]")
            println("Y range: [", minimum(my), ", ", maximum(my), "]") 
            println("Z range: [", minimum(mz_2d), ", ", maximum(mz_2d), "]")
        else
            rethrow(e)
        end
    end


    
    return nothing
end
