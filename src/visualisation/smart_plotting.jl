# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

if Threads.nthreads() > 1
    include("remote_plotting.jl") 
    init_plotting()  # This sets up workers and remote plotting capabilities   
end

"""
    smart_plot_flow_field(wf, X, Y, Z, vis; msr=1, plt=nothing)

High-level plotting function that automatically dispatches to either parallel or 
sequential plotting based on the number of available threads and processes.

# Arguments
- `wf`: WindFarm object
- `X`, `Y`, `Z`: Flow field coordinate arrays
- `vis`: Visualization settings
- `msr`: Measurement type (1=velocity reduction, 2=turbulence, 3=wind speed)
- `plt`: Matplotlib pyplot instance (only used in sequential mode)

# Returns
- Future object if using parallel execution, nothing otherwise
"""
function smart_plot_flow_field(wf, X, Y, Z, vis; msr=1, plt=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Use parallel plotting with remote worker
        return @spawnat 2 plot_flow_field(wf, X, Y, Z, vis; msr=msr)
    else
        # Use sequential plotting
        if plt === nothing
            error("plt argument is required for sequential plotting")
        end
        return plotFlowField(plt, wf, X, Y, Z, vis; msr=msr)
    end
end

"""
    smart_plot_measurements(wf, md, vis; separated=true, plt=nothing)

High-level measurements plotting function that automatically dispatches to either 
parallel or sequential plotting based on the number of available threads and processes.

# Arguments
- `wf`: WindFarm object
- `md`: Measurement data
- `vis`: Visualization settings
- `separated`: Whether to use separated subplots
- `plt`: Matplotlib pyplot instance (only used in sequential mode)

# Returns
- Future object if using parallel execution, nothing otherwise
"""
function smart_plot_measurements(wf, md, vis; separated=true, plt=nothing)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Use parallel plotting with remote worker
        return @spawnat 2 plot_measurements(wf, md, vis; separated=separated)
    else
        # Use sequential plotting
        if plt === nothing
            error("plt argument is required for sequential plotting")
        end
        return plotMeasurements(plt, wf, md, vis; separated=separated)
    end
end

