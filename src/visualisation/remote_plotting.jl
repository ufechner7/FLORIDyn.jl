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