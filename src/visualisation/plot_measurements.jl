# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    plotMeasurements(plt, wf::WindFarm, md::DataFrame, vis::Vis; separated=false) -> Nothing

Plot foreign reduction measurements from FLORIDyn simulation data.

# Arguments
- `plt`: Plotting package (e.g., ControlPlots)
- `wf::WindFarm`: Wind farm object with field `nT` (number of turbines). See [`WindFarm`](@ref)
- `md::DataFrame`: Measurements DataFrame containing time series data with columns:
  - `Time`: Time series data [s]
  - `ForeignReduction`: Foreign reduction percentage data [%]
- `vis::Vis`: Visualization settings including unit_test parameter. See [`Vis`](@ref)
- `separated::Bool`: Whether to use separated subplot layout (default: false)

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
plotMeasurements(plt, wind_farm, measurements_df, vis)

# Plot foreign reduction with separated subplots
plotMeasurements(plt, wind_farm, measurements_df, vis; separated=true)
```

# See Also
- [`plotFlowField`](@ref): For flow field visualization
- [`getMeasurements`](@ref): For generating measurement data
"""
function plotMeasurements(plt, wf::WindFarm, md::DataFrame, vis::Vis; separated=false)
    local fig

    # Subtract start time
    timeFDyn = md.Time .- md.Time[1]

    title="Foreign Reduction"
    size = 1
    
    # Foreign Reduction plotting
    if separated
        lay = get_layout(wf.nT)
        
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
    end
    # plt.show(fig)
    # Save plot to video or output folder if requested
    if vis.save && !vis.unit_test
        directory = vis.output_path
        
        # Generate filename with measurement type and time information
        msr = 1
        msr_name = msr == 1 ? "msr_velocity_reduction" : msr == 2 ? "msr_added_turbulence" : "msr_wind_speed"
        filename = joinpath(directory, "$(msr_name).png")
        
        # Save the current figure
        try
            plt.savefig(filename, dpi=150, bbox_inches="tight", pad_inches=0.1, facecolor="white")
            if !vis.unit_test
            end
        catch e
            @warn "Failed to save plot: $e"
        end
    end

    if vis.unit_test
        plt.pause(1.0)
        plt.close(fig)
    end
    plt.pause(0.01)
    return nothing
end

"""
    get_layout(nT::Int) -> (Int, Int)

Calculate optimal subplot layout (rows, columns) for a given number of turbines.

# Arguments
- `nT::Int`: Number of turbines

# Returns
- `Tuple{Int, Int}`: (rows, columns) for subplot arrangement

# Description
Determines the most square-like arrangement of subplots to accommodate `nT` plots.
"""
function get_layout(nT::Int)
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

