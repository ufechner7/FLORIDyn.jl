# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function init_plotting()
    # Only add a worker if we don't have any dedicated worker processes
    if nprocs() < 2  # nprocs() counts main + workers, so < 2 means no dedicated workers
        println("No dedicated workers found, adding 1 worker...")
        @time addprocs(1)
        @eval @everywhere using ControlPlots  # Ensure ControlPlots is available on all workers
        @eval @everywhere using FLORIDyn      # Ensure FLORIDyn (including WindFarm) is available on all workers
        @eval @everywhere using DataFrames    # Ensure DataFrames are available on all workers
         
    else
        println("Already have $(nprocs()-1) dedicated worker(s), total processes: $(nprocs())")
        println("Workers: $(workers())")
        
        # Ensure ControlPlots and FLORIDyn are loaded on existing workers
        @eval @everywhere using ControlPlots
        @eval @everywhere using FLORIDyn
    end
    nothing
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