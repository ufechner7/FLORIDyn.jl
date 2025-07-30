# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

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
            label = "Relative Wind Speed [%]"
            mz_2d .*= 100
            vmin = minimum(mz_2d); vmax = maximum(mz_2d);
        elseif msr == 2
            figure = "Added Turbulence"
            label = "Added Turbulence [%]"
            vmin = 0.0; vmax = maximum(mz_2d);
        elseif msr == 3
            figure = "Effective Wind Speed"
            vmin = 2.0; vmax = 10.0;
            label = L"Wind speed~[ms^{-1}]"
        end
        title = figure
        size = 0.84
        fig = plt.figure(figure, figsize=(7.25size, 6size))
        n=40
        levels = range(vmin, stop=vmax, length=n+1)
        plt.contourf(my, mx, mz_2d, n; levels, cmap="inferno") # 40 levels, no lines
        cb = plt.colorbar()
        cb.set_label(label, labelpad=3)

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

