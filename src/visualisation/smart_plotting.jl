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
    plot_rmt(X, Ys; xlabel="", ylabel="", labels=nothing, xlims=nothing, ylims=nothing, 
             ann=nothing, scatter=false, title="", fig="", ysize=14, pltctrl=nothing)

Smart plotting function for remote or sequential execution based on threading/processing environment.

This function automatically chooses between parallel plotting (using remote workers) and 
sequential plotting (using ControlPlots) based on the current Julia threading and 
multiprocessing configuration.

# Arguments
- `X`: X-axis data (typically a vector)
- `Ys`: Y-axis data. Can be a single vector or vector of vectors for multiple series
- `xlabel::String=""`: Label for the X-axis
- `ylabel::String=""`: Label for the Y-axis  
- `labels::Union{Nothing,Vector{String}}=nothing`: Labels for each data series (for legend)
- `xlims::Union{Nothing,Tuple}=nothing`: X-axis limits as (xmin, xmax)
- `ylims::Union{Nothing,Tuple}=nothing`: Y-axis limits as (ymin, ymax)
- `ann::Union{Nothing,Vector}=nothing`: Annotations to add to the plot
- `scatter::Bool=false`: Whether to create a scatter plot instead of line plot
- `title::String=""`: Plot title
- `fig::String=""`: Figure window title/identifier
- `ysize::Int=14`: Size of y-axis labels in points
- `pltctrl::Union{Nothing,Module}=nothing`: ControlPlots module for sequential plotting

# Behavior
- **Multithreaded + Multiprocessing**: Uses `@spawnat 2 Main.rmt_plot()` for parallel execution
- **Single-threaded or Single-process**: Uses `pltctrl.plot()` for sequential execution

# Threading Requirements
- **Parallel mode**: Requires `Threads.nthreads() > 1`, `nprocs() > 1`, and `pltctrl === nothing`
- **Sequential mode**: Requires `pltctrl` to be the ControlPlots module

# Examples
```julia
# Sequential plotting (single-threaded)
using ControlPlots
wind_dirs = [180.0, 200.0, 220.0, 240.0]
powers = [0.85, 0.92, 0.88, 0.76]
plot_rmt(wind_dirs, powers; xlabel="Wind Direction (deg)", ylabel="Relative Power", 
         title="Power vs Wind Direction", pltctrl=ControlPlots)

# Multiple data series
power_data = [[0.85, 0.92, 0.88, 0.76], [0.80, 0.89, 0.85, 0.72]]
plot_rmt(wind_dirs, power_data; xlabel="Wind Direction (deg)", ylabel="Relative Power",
         labels=["Configuration A", "Configuration B"], pltctrl=ControlPlots)

# Parallel plotting (multithreaded + multiprocessing)
# Automatically uses remote worker when threads > 1 and processes > 1
plot_rmt(wind_dirs, powers; xlabel="Wind Direction (deg)", ylabel="Relative Power")
```

# Notes
- The function returns `nothing` and displays the plot as a side effect
- In parallel mode, plotting occurs on worker process 2 via `@spawnat`
- In sequential mode, the plot is displayed using `display(p)` 
- Error handling ensures appropriate `pltctrl` parameter for sequential execution

# See Also
[`plot_x`](@ref), [`plot_flow_field`](@ref), [`plot_measurements`](@ref)
"""
function plot_rmt(X, Ys...; xlabel="", ylabel="", ylabels=nothing, labels=nothing, xlims=nothing, ylims=nothing, ann=nothing, 
    scatter=false, title="", fig="", ysize=14, pltctrl=nothing)

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