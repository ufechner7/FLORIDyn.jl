# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example of a simple Julia script to demonstrate the use of Distributed for plotting.

using Distributed

# Only add a worker if we don't have any dedicated worker processes
if nprocs() < 2  # nprocs() counts main + workers, so < 2 means no dedicated workers
    println("No dedicated workers found, adding 1 worker...")
    @time addprocs(1)
else
    println("Already have $(nprocs()-1) dedicated worker(s), total processes: $(nprocs())")
    println("Workers: $(workers())")
end

@everywhere using ControlPlots  # Ensure ControlPlots is available on all workers
using ControlPlots

@spawnat 2 display(plot(rand(3)))

nothing
