# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example of a simple Julia script to demonstrate the use of DistributedNext for plotting.

using DistributedNext, Timers, ControlPlots, FLORIDyn

tic()
include("../src/visualisation/remote_plotting.jl") 
init_plotting()
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
