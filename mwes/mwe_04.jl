# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example of a simple Julia script to demonstrate the use of Distributed for plotting.

using Distributed, Timers

# Only add a worker if we don't have any dedicated worker processes
if nprocs() < 2  # nprocs() counts main + workers, so < 2 means no dedicated workers
    println("No dedicated workers found, adding 1 worker...")
    @time addprocs(1)
else
    println("Already have $(nprocs()-1) dedicated worker(s), total processes: $(nprocs())")
    println("Workers: $(workers())")
end

tic()
@everywhere using ControlPlots  # Ensure ControlPlots is available on all workers
@everywhere using FLORIDyn      # Ensure FLORIDyn (including WindFarm) is available on all workers
using ControlPlots
using FLORIDyn
toc()

# Test that WindFarm is available on all workers
println("Testing WindFarm availability on workers...")
for worker_id in workers()
    try
        # Test if WindFarm type is available on the worker
        result = @spawnat worker_id typeof(FLORIDyn.WindFarm)
        println("Worker $worker_id: WindFarm type = ", fetch(result))
    catch e
        println("Worker $worker_id: ERROR - ", e)
    end
end

X = 0:0.1:10
Y = sin.(X)

REMOTE = true  # Set to true if you want to plot on a remote worker
if REMOTE
    # If REMOTE, plot on a specific worker (e.g., worker 2)
    @time @spawnat 2 display(plot(X, Y, fig="Sine Function Plot", xlabel="X-axis", ylabel="Y-axis"))
else
    # If not REMOTE, just plot on the main process
    @time display(plot(X, Y, fig="Sine Function Plot", xlabel="X-axis", ylabel="Y-axis"))
end

# Results
# remote:   0.000199 seconds (185 allocations: 12.336 KiB)
# main:     0.077712 seconds (365 allocations: 10.672 KiB)

nothing
