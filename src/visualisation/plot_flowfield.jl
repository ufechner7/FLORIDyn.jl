# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    PlotState

Mutable struct to hold animation state for flow field plotting.

# Fields
- `fig`: Matplotlib figure object
- `ax`: Matplotlib axes object  
- `cb`: Colorbar object
- `contour_collection`: Collection from contourf for updating data
- `turbine_lines`: Vector of line objects for turbine rotors
- `op_scatter1`: Scatter plot object for all operational points
- `op_scatter2`: Scatter plot object for every 10th operational point
- `title_obj`: Title text object for updating
- `figure_name`: Name of the figure window
- `label`: Colorbar label text
- `lev_min`: Minimum value for color scale
- `lev_max`: Maximum value for color scale
- `levels`: Contour levels
"""
mutable struct PlotState
    fig
    ax
    cb
    contour_collection
    turbine_lines::Vector
    op_scatter1
    op_scatter2
    title_obj
    figure_name::String
    label::String
    lev_min::Float64
    lev_max::Float64
    levels
end

"""
    plotFlowField(state::Union{Nothing, PlotState}, plt, wf, mx, my, mz; msr=3, unit_test=false)

Plot a 2D contour of the flow field data with support for animation.

# Arguments
- `state::Union{Nothing, PlotState}`: Animation state object. Pass `nothing` for the first call, 
  then pass the returned state for subsequent calls to maintain the same figure and layout.
- `plt`: Plotting package (e.g., ControlPlots)
- `wf`: Wind farm object containing turbine data
- `mx::Matrix`: X-coordinate grid
- `my::Matrix`: Y-coordinate grid  
- `mz::Array{Float64,3}`: 3D array of measurements with dimensions (rows, cols, nM)
- `msr::Int`: Which measurement to plot (1, 2, or 3). Default is 3.
- `unit_test::Bool`: Whether to automatically close plots for testing. Default is false.

# Returns
- `state::PlotState`: Updated or newly created plot state for use in subsequent calls

# Description
This function supports creating animations by maintaining plot state across multiple calls:

## First Call (state = nothing)
- Creates new figure, axes, colorbar, and all plot elements
- Initializes and returns a PlotState object

## Subsequent Calls (state = PlotState)
- Updates existing contour data, turbine positions, and operational points
- Reuses the same figure and layout for smooth animation

# Animation Example
```julia
using ControlPlots

# Initialize state (first frame)
state = nothing
for t in time_steps
    Z, X, Y = calcFlowField(settings, wind_farm, wind, floris)
    state = plotFlowField(state, plt, wind_farm, X, Y, Z; msr=1)
    plt.pause(0.1)  # Small delay for animation
end
```

# Measurement Types
- `msr=1`: Velocity reduction [%]
- `msr=2`: Added turbulence [%]  
- `msr=3`: Effective wind speed [m/s]

# Notes
- The function automatically handles coordinate system transformations for turbine orientations
- Operational points are displayed as white scatter points for reference
- Color scales are kept consistent across animation frames when using the same measurement type
- This function requires a plotting package like ControlPlots.jl to be loaded and available as `plt`
"""
function plotFlowField(state::Union{Nothing, PlotState}, plt, wf, mx, my, mz; msr=3, unit_test=false)
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
            figure_name = "Velocity Reduction"
            label = "Relative Wind Speed [%]"
            mz_2d .*= 100
            lev_min = min(vis.rel_v_min, minimum(mz_2d)); lev_max = maximum(mz_2d);
        elseif msr == 2
            figure_name = "Added Turbulence"
            label = "Added Turbulence [%]"
            lev_min = 0.0; lev_max = maximum(mz_2d);
        elseif msr == 3
            figure_name = "Effective Wind Speed"
            lev_min = 2.0; lev_max = 10.0;
            label = L"Wind speed~[ms^{-1}]"
        end
        title = figure_name
        
        # Initialize or update plot state
        if state === nothing
            # First call - create new figure and all elements
            size = 0.84
            fig = plt.figure(figure_name, figsize=(7.25size, 6size))
            ax = plt.gca()
            n = 40
            levels = range(lev_min, stop=lev_max, length=n+1)
            contour_collection = plt.contourf(my, mx, mz_2d, n; levels, cmap="inferno")
            cb = plt.colorbar()
            cb.set_label(label, labelpad=3)
            
            # Initialize empty containers for dynamic elements
            turbine_lines = []
            op_scatter1 = nothing
            op_scatter2 = nothing
            title_obj = plt.title(title)
            
            # Create state object
            state = PlotState(fig, ax, cb, contour_collection, turbine_lines, op_scatter1, op_scatter2, 
                             title_obj, figure_name, label, lev_min, lev_max, levels)
        else
            # Subsequent calls - update existing plot
            plt.figure(state.figure_name)
            
            # Clear previous dynamic elements
            for line in state.turbine_lines
                try
                    line.remove()
                catch
                    # Line might already be removed
                end
            end
            empty!(state.turbine_lines)
            
            if state.op_scatter1 !== nothing
                try
                    state.op_scatter1.remove()
                catch
                end
            end
            if state.op_scatter2 !== nothing
                try
                    state.op_scatter2.remove()
                catch
                end
            end
            
            # Update contour data
            for collection in state.contour_collection.collections
                collection.remove()
            end
            
            # Create new contour with updated data
            state.contour_collection = plt.contourf(my, mx, mz_2d, 40; levels=state.levels, cmap="inferno")
        end

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

            # Plot turbine line and store reference
            ax = plt.gca()
            line = ax.plot(rot_pos[1, :], rot_pos[2, :], [20, 20], color="k", linewidth=3)[1]
            push!(state.turbine_lines, line)
        end
        
        # Plot the OPs
        # Plot all points with size 2 and white filled marker
        state.op_scatter1 = plt.scatter(wf.States_OP[:, 1], wf.States_OP[:, 2], s=2, color="white", marker="o")
        
        # Plot every 10th point with size 6 and white filled marker
        state.op_scatter2 = plt.scatter(wf.States_OP[1:10:end, 1], wf.States_OP[1:10:end, 2], s=6, color="white", marker="o")
        
        plt.xlim(minimum(mx), maximum(mx))
        plt.ylim(minimum(mx), maximum(mx))
        state.title_obj.set_text(title)
        plt.xlabel("West-East [m]")
        plt.ylabel("South-North [m]")
        plt.tight_layout()
        
        # Force display update for animation
        plt.draw()
        
        if !unit_test
            println("Contour plot updated successfully")
        end
        
        if unit_test
            plt.pause(2)
            plt.close(state.fig)
            return nothing
        end
    catch e
        if isa(e, UndefVarError) && e.var == :plt
            @warn "ControlPlots not available. Please load ControlPlots.jl first with: using ControlPlots"
            println("Data shape: ", size(mz_2d))
            println("X range: [", minimum(mx), ", ", maximum(mx), "]")
            println("Y range: [", minimum(my), ", ", maximum(my), "]") 
            println("Z range: [", minimum(mz_2d), ", ", maximum(mz_2d), "]")
            return state
        else
            rethrow(e)
        end
    end
    return state
end

"""
    plotFlowField(plt, wf, mx, my, mz; msr=3, unit_test=false)

Compatibility method for the original plotFlowField interface.

This method provides backward compatibility by calling the new state-based version 
with `state=nothing`, effectively creating a single plot without animation support.

# Arguments
- `plt`: Plotting package (e.g., ControlPlots)
- `wf`: Wind farm object containing turbine data
- `mx::Matrix`: X-coordinate grid
- `my::Matrix`: Y-coordinate grid  
- `mz::Array{Float64,3}`: 3D array of measurements with dimensions (rows, cols, nM)
- `msr::Int`: Which measurement to plot (1, 2, or 3). Default is 3.
- `unit_test::Bool`: Whether to automatically close plots for testing. Default is false.

# Returns
- `nothing`: For compatibility with the original interface

# Note
This method is provided for backward compatibility. For animation support, 
use the new interface with explicit state management.
"""
function plotFlowField(plt, wf, mx, my, mz; msr=3, unit_test=false)
    plotFlowField(nothing, plt, wf, mx, my, mz; msr=msr, unit_test=unit_test)
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
            fig.subplots_adjust(wspace=0.295)
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

