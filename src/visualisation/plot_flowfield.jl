# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

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
    
    # Single-threaded loop (can be parallelized with @threads or Distributed.@distributed)
    for iGP in 1:length(mx)
        xGP = mx[iGP]
        yGP = my[iGP]
        
        # Update only the grid point position in the pre-allocated arrays
        GP.posBase[end, :] = [xGP, yGP, 0.0]  # Update grid point base position
        GP.posNac[end, :] = [0.0, 0.0, zh]    # Update grid point nacelle position
        
        # Reset the grid point state (other turbine states remain unchanged)
        GP.States_T[end, :] .= 0.0
        
        # Recalculate interpolated OPs for the updated geometry
        GP.intOPs = interpolateOPs(GP)

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

This function requires a plotting package like ControlPlots.jl to be loaded and available as `plt`.
"""
function plotFlowField(plt, wf, mx, my, mz; msr=3, unit_test=false)
    # Extract the 2D slice for the specified measurement
    if msr > size(mz, 3)
        error("msr ($msr) exceeds number of measurements ($(size(mz, 3)))")
    end
    
    # Get the 2D slice
    mz_2d = mz[:, :, msr]
    
    # Try to use ControlPlots if available
    try
        # This will work if ControlPlots is loaded and plt is available
        if msr == 1
            figure = "Velocity Reduction"
        elseif msr == 2
            figure = "Added Turbulence"
        elseif msr == 3
            figure = "Effective Wind Speed"
        end
        title = figure
        size = 0.84
        fig = plt.figure(figure, figsize=(7.25size, 6size))
        vmin = 2.0; vmax = 10.0; n=40
        levels = range(vmin, stop=vmax, length=n+1)
        plt.contourf(my, mx, mz_2d, n; levels, cmap="inferno") # 40 levels, no lines
        cb = plt.colorbar()
        cb.set_label(L"Wind speed~[ms^{-1}]", labelpad=3)

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
        # Plot all points with size 2 and white filled marker
        n = 20
        plt.scatter(wf.States_OP[:, 1], wf.States_OP[:, 2], s=2, color="white", marker="o")
        plt.xlim(minimum(mx), maximum(mx))
        plt.ylim(minimum(mx), maximum(mx))

        # Plot every 10th point with size 6 and white filled marker
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
            @warn "ControlPlots not available. Please load ControlPlots.jl first with: using ControlPlots"
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

"""
    plotMeasurements(plt, wf::WindFarm, md::DataFrame; separated=false, unit_test=false) -> Nothing

Plot foreign reduction measurements from FLORIDyn simulation data.

# Arguments
- `plt`: Plotting package (e.g., ControlPlots)
- `wf::WindFarm`: Wind farm object with field `nT` (number of turbines)
- `md::DataFrame`: Measurements DataFrame containing time series data with columns:
  - `Time`: Time series data [s]
  - `ForeignReduction`: Foreign reduction percentage data [%]
- `separated::Bool`: Whether to use separated subplot layout (default: false)
- `unit_test::Bool`: Whether to close plots automatically for testing (default: false)

# Returns
- `nothing`

# Description
This function creates time series plots of foreign reduction measurements from FLORIDyn simulations. It handles:
1. Time normalization by subtracting the start time
2. Foreign reduction plots in either separated (subplot) or combined layout

# Plotting Modes
- **Separated mode** (`separated=true`): Creates individual subplots for each turbine
- **Combined mode** (`separated=false`): Plots all turbines on a single figure with different colors

# Example
```julia
using ControlPlots

# Plot foreign reduction for all turbines in combined mode
plotMeasurements(plt, wind_farm, measurements_df)

# Plot foreign reduction with separated subplots
plotMeasurements(plt, wind_farm, measurements_df; separated=true)
```

# See Also
- [`plotFlowField`](@ref): For flow field visualization
- [`getMeasurements`](@ref): For generating measurement data
"""
function plotMeasurements(plt, wf::WindFarm, md::DataFrame; separated=false, unit_test=false)
    local fig

    # Subtract start time
    timeFDyn = md.Time .- md.Time[1]

    title="Foreign Reduction"
    size = 1
    
    # Foreign Reduction plotting
    if separated
        lay = getLayout(wf.nT)
        
        # Calculate y-axis limits
        foreign_red_data = md[!, "ForeignReduction"]
        y_lim = [minimum(foreign_red_data), maximum(foreign_red_data)]
        y_range = y_lim[2] - y_lim[1]
        y_lim = y_lim .+ [-0.1, 0.1] * max(y_range, 0.5)
        
        fig = plt.figure(title, figsize=(10size, 6size))
        c = plt.get_cmap("inferno")(0.5)  # Single color for separated plots
        
        for iT in 1:wf.nT
            plt.subplot(lay[1], lay[2], iT)
            plt.plot(
                timeFDyn[iT:wf.nT:end],
                foreign_red_data[iT:wf.nT:end],
                linewidth=2, color=c
            )
            plt.grid(true)
            plt.title("Turbine $(iT)")
            plt.xlim(max(timeFDyn[1], 0), timeFDyn[end])
            plt.ylim(y_lim...)
            plt.xlabel("Time [s]")
            plt.ylabel("Foreign Reduction [%]")
            plt.tight_layout()   
        end
        println("Foreign reduction plot created successfully")
    else
        fig = plt.figure(title*" - Line Plot")
        colors = plt.get_cmap("inferno")(range(0, 1, length=wf.nT))
        
        for iT in 1:wf.nT
            plt.plot(
                timeFDyn[iT:wf.nT:end],
                md.ForeignReduction[iT:wf.nT:end],
                linewidth=2, color=colors[iT, :]
            )
        end
        plt.grid(true)
        plt.title("Foreign Reduction [%]")
        plt.xlim(max(timeFDyn[1], 0), timeFDyn[end])
        plt.xlabel("Time [s]")
        plt.ylabel("Foreign Reduction [%]")
        println("Foreign reduction line plot created successfully")
    end
    if unit_test
        plt.pause(2)
        plt.close(fig)
    end
        
    return nothing
end

"""
    getLayout(nT::Int) -> Tuple{Int, Int}

Calculate optimal subplot layout (rows, columns) for a given number of turbines.

# Arguments
- `nT::Int`: Number of turbines

# Returns
- `Tuple{Int, Int}`: (rows, columns) for subplot arrangement

# Description
Determines the most square-like arrangement of subplots to accommodate `nT` plots.
"""
function getLayout(nT::Int)
    if nT <= 0
        return (1, 1)
    elseif nT == 1
        return (1, 1)
    elseif nT <= 4
        return (2, 2)
    elseif nT <= 6
        return (2, 3)
    elseif nT <= 9
        return (3, 3)
    elseif nT <= 12
        return (3, 4)
    elseif nT <= 16
        return (4, 4)
    else
        # For larger numbers, calculate a roughly square layout
        cols = ceil(Int, sqrt(nT))
        rows = ceil(Int, nT / cols)
        return (rows, cols)
    end
end

