# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    plot_flow_field(wf, X, Y, Z, vis; msr=1, plt=nothing) -> Nothing

High-level plotting function that automatically dispatches to either parallel or 
sequential plotting based on the number of available threads and processes.

# Arguments
- `wf`: WindFarm object
- `X`, `Y`, `Z`: Flow field coordinate arrays
- `vis`: Visualization settings
- `msr`: Measurement type (1=velocity reduction, 2=turbulence, 3=wind speed)
- `plt`: Matplotlib PyPlot instance (only used in sequential mode)

# Returns
- nothing
"""
function plot_flow_field(wf, X, Y, Z, vis; msr=1, plt=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Use parallel plotting with remote worker
        @spawnat 2 Main.rmt_plot_flow_field(wf, X, Y, Z, vis; msr=msr)
    else
        # Use sequential plotting
        if plt === nothing
            error("plt argument is required for sequential plotting")
        end
        plotFlowField(plt, wf, X, Y, Z, vis; msr=msr)
    end
    nothing
end

"""
    plot_measurements(wf, md, vis; separated=true, plt=nothing) -> Nothing

High-level measurements plotting function that automatically dispatches to either 
parallel or sequential plotting based on the number of available threads and processes.

# Arguments
- `wf`: WindFarm object
- `md`: Measurement data
- `vis`: Visualization settings
- `separated`: Whether to use separated subplots
- `plt`: Matplotlib PyPlot instance (only used in sequential mode)

# Returns
- nothing

# See Also
- [`plotMeasurements`](@ref): The underlying plotting function used in sequential mode
"""
function plot_measurements(wf, md, vis; separated=true, plt=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Use parallel plotting with remote worker
        @spawnat 2 Main.rmt_plot_measurements(wf, md, vis; separated=separated)
    else
        # Use sequential plotting
        if plt === nothing
            error("plt argument is required for sequential plotting")
        end
        plotMeasurements(plt, wf, md, vis; separated=separated)        
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