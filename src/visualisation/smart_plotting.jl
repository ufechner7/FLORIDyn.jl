# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt=nothing, 
                    fig=nothing) -> Nothing

High-level plotting function that automatically dispatches to either parallel or 
sequential plotting based on the number of available threads and processes.

# Arguments
- `wf`: WindFarm object
- `X`, `Y`, `Z`: Flow field coordinate arrays
- `vis`: Visualization settings
- `msr`: Measurement type, see: [MSR](@ref)
- `plt`: Matplotlib PyPlot instance (only used in sequential mode)
- `fig`: Figure name (optional)

# Returns
- nothing
"""
function plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt=nothing, fig=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Use parallel plotting with remote worker
        @spawnat 2 Main.rmt_plot_flow_field(wf, X, Y, Z, vis; msr, fig)
    else
        # Use sequential plotting
        if plt === nothing
            error("plt argument is required for sequential plotting")
        end
        plotFlowField(plt, wf, X, Y, Z, vis; msr, fig)
    end
    nothing
end

"""
    plot_measurements(wf, md, vis; separated=true, msr=VelReduction, plt=nothing) -> Nothing

High-level measurements plotting function that automatically dispatches to either 
parallel or sequential plotting based on the number of available threads and processes.

# Arguments
- `wf`: WindFarm object
- `md`: Measurement data
- `vis`: Visualization settings
- `separated`: Whether to use separated subplots
- `msr`: Measurement type, see: [MSR](@ref) 
- `plt`: Matplotlib PyPlot instance (only used in sequential mode)

# Returns
- nothing

# See Also
- [`plotMeasurements`](@ref): The underlying plotting function used in sequential mode
"""
function plot_measurements(wf, md, vis; separated=true, msr::MSR=VelReduction, plt=nothing, pltctrl=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Parallel/remote plotting path (doesn't need local ControlPlots instance)
        @spawnat 2 Main.rmt_plot_measurements(wf, md, vis; separated, msr)
    else
        # Sequential plotting path with backward compatibility
        # Accept either the new ControlPlots module via pltctrl OR a legacy PyPlot-like plt.
        local_plt = nothing
        if pltctrl !== nothing
            # New style: module ControlPlots passed -> use its internal PyPlot handle
            try
                local_plt = pltctrl.plt
            catch
                error("pltctrl provided but no .plt field found; pass the ControlPlots module itself")
            end
        elseif plt !== nothing
            # Legacy style: a PyPlot-like object was passed directly
            local_plt = plt
        else
            error("Either pltctrl (ControlPlots module) or plt (legacy PyPlot handle) must be provided for sequential plotting")
        end
        plotMeasurements(local_plt, wf, md, vis; separated, msr)
    end
    return nothing
end

"""
    plot_x(times, plot_data...; ylabels=nothing, labels=nothing, fig="Wind Direction", 
           xlabel="rel_time [s]", ysize=10, bottom=0.02, legend_size=nothing, pltctrl=nothing, 
           loc=nothing) -> Nothing

High-level time series plotting function that automatically dispatches to either 
parallel or sequential plotting based on the number of available threads and processes.

# Arguments
- `times`:        Time vector for x-axis
- `plot_data...`: Variable number of data arrays to plot
- `ylabels`:      Labels for y-axes (optional)
- `labels`:       Labels for subplots (optional)
- `fig`:          Figure title (default: "Wind Direction")
- `xlabel`:       X-axis label (default: "rel_time [s]")
- `ysize`:        Size of the Y-axis labels in points (default: 10)
- `bottom`:       Bottom margin (default: 0.02)
- `legend_size`:  Legend font size in points (optional)
- `pltctrl`:      ControlPlots instance (only used in sequential mode)

# Returns
- nothing

# Description
When running with multiple threads and processes, it uses remote plotting 
capabilities via `rmt_plotx` ONLY if no local `pltctrl` is provided. If a `pltctrl`
argument is supplied (e.g. during unit tests with a mock), it forces the sequential
path so tests can observe side-effects without needing remote worker setup.

# Example
```julia
plot_x(times, data1, data2; ylabels=["Turbine 1", "Turbine 2"], 
       labels=["Wind Speed", "Power"], legend_size=8, pltctrl=pltctrl)
```

# See Also
- `plotx`: The underlying plotting function from ControlPlots used in sequential mode
"""
function plot_x(times, plot_data...; ylabels=nothing, labels=nothing, 
                fig="Wind Direction", xlabel="rel_time [s]", ysize=10, bottom=0.02, 
                legend_size=nothing, pltctrl=nothing, loc=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1 && pltctrl === nothing
        # Use parallel plotting with remote worker
        @spawnat 2 Main.rmt_plotx(times, plot_data...; ylabels=ylabels, labels=labels,
                                  fig=fig, xlabel=xlabel, ysize=ysize, bottom=bottom, legend_size=legend_size, loc=loc)
    else
        # Use sequential plotting
        if pltctrl === nothing
            error("pltctrl argument is required for sequential plotting (threads=$(Threads.nthreads()), procs=$(nprocs())). Pass the ControlPlots module as pltctrl keyword.")
        end
        p = pltctrl.plotx(times, plot_data...; ylabels=ylabels, labels=labels,
                           fig=fig, xlabel=xlabel, ysize=ysize, bottom=bottom, legend_size=legend_size, loc=loc)
        display(p)  # Ensure the plot is displayed
    end
    nothing
end

"""
    plot_rmt(X, Ys...; xlabel="", ylabel="", ylabels=nothing, labels=nothing, 
             xlims=nothing, ylims=nothing, ann=nothing, scatter=false, title="", 
             fig="", ysize=14, pltctrl=nothing) -> Nothing

Smart plotting function that automatically dispatches between remote parallel plotting 
and sequential plotting based on Julia's threading and multiprocessing environment.

This function provides a unified interface for creating line plots or scatter plots
that adapts to the available computational resources. It supports both single and 
multiple data series with comprehensive labeling and styling options.

# Arguments
- `X`: X-axis data vector (e.g., time points, wind directions, positions)
- `Ys...`: Variable number of Y-axis data vectors. Each should have the same length as `X`
- `xlabel::String=""`: Label for the X-axis
- `ylabel::String=""`: Label for the Y-axis (only valid for single Y axis)
- `ylabels::Union{Nothing,Vector{String}}=nothing`: Y-axis labels for two Y axis
- `labels::Union{Nothing,Vector{String}}=nothing`: Legend labels for each data series
- `xlims::Union{Nothing,Tuple{Real,Real}}=nothing`: X-axis limits as `(xmin, xmax)`
- `ylims::Union{Nothing,Tuple{Real,Real}}=nothing`: Y-axis limits as `(ymin, ymax)`
- `ann::Union{Nothing,Vector}=nothing`: Annotations to add to the plot
- `scatter::Bool=false`: Create scatter plot if `true`, line plot if `false`
- `title::String=""`: Plot title displayed at the top
- `fig::String=""`: Figure window title/identifier for window management
- `ysize::Int=14`: Font size for y-axis labels in points
- `pltctrl::Union{Nothing,Module}=nothing`: ControlPlots module instance for sequential plotting

# Returns
- `Nothing`: Function executes for side effects (displays plot)

# Execution Modes

## Parallel Mode (Remote Plotting)
**Conditions**: `Threads.nthreads() > 1` AND `nprocs() > 1` AND `pltctrl === nothing`
- Executes plotting on remote worker process 2 using `@spawnat 2 Main.rmt_plot(...)`
- Reduces main process computational load for large datasets

## Sequential Mode (Local Plotting)  
**Conditions**: Single-threaded OR single-process OR `pltctrl` provided
- Uses provided `pltctrl` (ControlPlots module) for local execution
- Required for unit testing and environments without multiprocessing
- Calls `pltctrl.plot(...)` with appropriate parameter combinations

# Parameter Validation
- **Length consistency**: All `Ys` vectors must match `X` length
- **Label compatibility**: `ylabel` cannot be used with multiple data series (use `ylabels`)
- **Required arguments**: `pltctrl` is mandatory for sequential mode

# Examples

## Single Data Series
```julia
using ControlPlots
time = 0:0.1:10
power = sin.(time) .+ 1
plot_rmt(time, power; 
         xlabel="Time (s)", ylabel="Power (MW)", 
         title="Turbine Power Output", pltctrl=ControlPlots)
```

## Multiple Data Series with Legend
```julia
wind_dirs = [180, 200, 220, 240, 260, 280]
power_t1 = [0.85, 0.92, 0.88, 0.76, 0.65, 0.58]
power_t2 = [0.80, 0.89, 0.85, 0.72, 0.61, 0.55]

plot_rmt(wind_dirs, power_t1, power_t2;
         xlabel="Wind Direction (°)", ylabels=["Power T1", "Power T2"],
         labels=["Turbine 1", "Turbine 2"], 
         title="Power vs Wind Direction",
         xlims=(175, 285), pltctrl=ControlPlots)
```

## Scatter Plot with Annotations
```julia
x_pos = [0, 500, 1000, 1500]
y_pos = [0, 200, -100, 300]
plot_rmt(x_pos, y_pos; 
         xlabel="X Position (m)", ylabel="Y Position (m)",
         scatter=true, title="Turbine Layout",
         ann=["T1", "T2", "T3", "T4"], pltctrl=ControlPlots)
```

## Automatic Parallel Execution
```julia
# When Julia started with: julia -p 2 -t 4
# Automatically uses remote plotting (no pltctrl needed)
plot_rmt(wind_speeds, power_curve; 
         xlabel="Wind Speed (m/s)", ylabel="Power (MW)",
         title="Power Curve")
```

# Error Handling
- `ArgumentError`: Length mismatch between `X` and any `Ys` vector
- `ArgumentError`: Using `ylabel` with multiple data series  
- `ErrorException`: Missing `pltctrl` in sequential mode

# Implementation Details
The function uses conditional dispatch based on Julia's threading and multiprocessing state:
1. Checks `Threads.nthreads() > 1 && nprocs() > 1 && pltctrl === nothing`
2. **If true**: Parallel execution via `@spawnat 2 Main.rmt_plot(...)`
3. **If false**: Sequential execution via `pltctrl.plot(...)` with parameter branching

Parameter branching in sequential mode:
- Single series: Uses `ylabel` parameter
- Multiple series without `ylabels`: Uses `ylabel` and `labels`  
- Multiple series with `ylabels`: Uses `ylabels`, `labels`, omits `ylabel`

# See Also
- [`plot_x`](@ref): Time series plotting with automatic subplot generation
- [`plot_flow_field`](@ref): 2D/3D flow field visualization
- [`plot_measurements`](@ref): Wind farm measurement data plotting
- [`ControlPlots.plot`]: Underlying plotting function for sequential mode
"""
function plot_rmt(X, Ys...; xlabel="", ylabel="", ylabels=nothing, labels=nothing, xlims=nothing, ylims=nothing, ann=nothing, 
    scatter=false, title="", fig="", ysize=14, pltctrl=nothing)
    
    # Parameter validation: Ensure X and each Y in Ys have compatible dimensions
    for (i, Y) in pairs(Ys...)
        if length(X) != length(Y)
            throw(ArgumentError("Length of X ($(length(X))) does not match length of Ys[$i] ($(length(Y)))."))
        end
    end
    # Validate ylabel vs multiple Y series
    if ylabel != "" && length(Ys) > 1
        throw(ArgumentError("Cannot use ylabel with multiple Y series (detected $(length(Ys))). Use ylabels instead, e.g. ylabels=[\"Series 1\", \"Series 2\"]."))
    end

    if Threads.nthreads() > 1 && nprocs() > 1 && pltctrl === nothing
        # Use parallel plotting with remote worker

        @spawnat 2 Main.rmt_plot(X, Ys...; xlabel, ylabel, ylabels, labels, xlims, ylims, ann, scatter, title, fig, ysize)
    else
        # Use sequential plotting
        if pltctrl === nothing
            error("pltctrl argument is required for sequential plotting (threads=$(Threads.nthreads()), procs=$(nprocs())). Pass the ControlPlots module as pltctrl keyword.")
        end
        if isnothing(labels) && length(Ys) == 1
            p = pltctrl.plot(X, Ys...; xlabel, ylabel, xlims, ylims, ann, scatter, title, fig, ysize)
        elseif isnothing(ylabels)
            p = pltctrl.plot(X, Ys...; xlabel, ylabel, labels, xlims, ylims, ann, scatter, title, fig, ysize)
        else
            p = pltctrl.plot(X, Ys...; xlabel, ylabels, labels, title, fig, ysize)
        end
        
        display(p)  # Ensure the plot is displayed
    end
    nothing
end

"""
    close_all(plt)

Close all matplotlib figure windows.

This function automatically dispatches to either parallel or sequential plotting
based on the number of available threads and processes.

# Arguments
- `plt`: Matplotlib PyPlot instance (only used in sequential mode)

# Description
When running with multiple threads and processes, it uses remote plotting 
capabilities to close all figures on the remote worker. Otherwise, it directly
calls `plt.close("all")` to close all figures in the current process.
"""
function close_all(plt)
    if Threads.nthreads() > 1 && nprocs() > 1
        @spawnat 2 Main.rmt_close_all()
    else
        plt.close("all")
    end
end