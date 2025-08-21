# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt=nothing, fig=nothing) -> Nothing

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
function plot_measurements(wf, md, vis; separated=true, msr::MSR=VelReduction, plt=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Use parallel plotting with remote worker
        @spawnat 2 Main.rmt_plot_measurements(wf, md, vis; separated=separated, msr)
    else
        # Use sequential plotting
        if plt === nothing
            error("plt argument is required for sequential plotting")
        end
        plotMeasurements(plt, wf, md, vis; separated=separated, msr)        
    end
    nothing
end

"""
    plot_x(times, plot_data...; ylabels=nothing, labels=nothing, fig="Wind Direction", 
           xlabel="rel_time [s]", ysize=10, bottom=0.02, plt=nothing) -> Nothing

High-level time series plotting function that automatically dispatches to either 
parallel or sequential plotting based on the number of available threads and processes.

# Arguments
- `times`: Time vector for x-axis
- `plot_data...`: Variable number of data arrays to plot
- `ylabels`: Labels for y-axes (optional)
- `labels`: Labels for subplots (optional)
- `fig`: Figure title (default: "Wind Direction")
- `xlabel`: X-axis label (default: "rel_time [s]")
- `ysize`: Figure height (default: 10)
- `bottom`: Bottom margin (default: 0.02)
- `plt`: Matplotlib PyPlot instance (only used in sequential mode)

# Returns
- nothing

# Description
When running with multiple threads and processes, it uses remote plotting 
capabilities via `rmt_plotx`. Otherwise, it directly calls `plotx` with the 
provided pyplot instance.

# Example
```julia
plot_x(times, data1, data2; ylabels=["Turbine 1", "Turbine 2"], 
       labels=["Wind Speed", "Power"], plt=plt)
```

# See Also
- [`plotx`](@ref): The underlying plotting function used in sequential mode
"""
function plot_x(times, plot_data...; ylabels=nothing, labels=nothing, 
                fig="Wind Direction", xlabel="rel_time [s]", ysize=10, bottom=0.02, plt=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Use parallel plotting with remote worker
        @spawnat 2 Main.rmt_plotx(times, plot_data...; ylabels=ylabels, labels=labels,
                                  fig=fig, xlabel=xlabel, ysize=ysize, bottom=bottom)
    else
        # Use sequential plotting
        if plt === nothing
            error("plt argument is required for sequential plotting")
        end
        p=plt.plotx(times, plot_data...; ylabels=ylabels, labels=labels,
              fig=fig, xlabel=xlabel, ysize=ysize, bottom=bottom)
        display(p)  # Ensure the plot is displayed in interactive mode
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