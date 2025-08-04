# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

if Threads.nthreads() > 1
    include("remote_plotting.jl") 
    init_plotting()  # This sets up workers and remote plotting capabilities   
end

"""
    plot_measurements(wf, md, vis; separated=true, plt=nothing)

High-level measurements plotting function that automatically dispatches to either 
parallel or sequential plotting based on the number of available threads and processes.

# Arguments
- `wf`: WindFarm object
- `md`: Measurement data
- `vis`: Visualization settings
- `separated`: Whether to use separated subplots
- `plt`: Matplotlib PyPlot instance (only used in sequential mode)

# Returns
- Future object if using parallel execution, nothing otherwise
"""
function plot_measurements(wf, md, vis; separated=true, plt=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Use parallel plotting with remote worker
        return @spawnat 2 rmt_plot_measurements(wf, md, vis; separated=separated)
    else
        # Use sequential plotting
        if plt === nothing
            error("plt argument is required for sequential plotting")
        end
        return plotMeasurements(plt, wf, md, vis; separated=separated)
    end
end
