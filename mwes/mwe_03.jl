# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example of a simple Julia script to demonstrate the use of DistributedNext for parallel processing
# This script initializes a worker, defines a function, and runs it on the worker

using Distributed

# Only add a worker if we don't have any dedicated worker processes
if nprocs() < 2  # nprocs() counts main + workers, so < 2 means no dedicated workers
    println("No dedicated workers found, adding 1 worker...")
    @time addprocs(1)
else
    println("Already have $(nprocs()-1) dedicated worker(s), total processes: $(nprocs())")
    println("Workers: $(workers())")
end

# Define function on all processes (main + workers)
@everywhere function identity_func(x)
    return x
end

function pltr(data)
    worker_id = workers()[1]
    remote_do(identity_func, worker_id, data)
end

data = rand(UInt8, 1024^2)  # 1MB random data
worker_id = workers()[1]  # Get the first available worker ID
println("Using worker $worker_id")

# Pre-compile by calling once
pltr(data)
@time pltr(data)
nothing
