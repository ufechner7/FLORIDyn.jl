# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, BenchmarkTools, Parameters

import Base: show

settings_file = "data/2021_9T_Data.yaml"
wind, sim, con, floris, floridyn, ta = setup(settings_file)
set = Settings(wind, sim, con)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf = initSimulation(wf, sim)

# Initialize dependencies for all turbines
wf.dep = [
    Int64[4, 7],
    Int64[5, 7, 8],
    Int64[6, 8, 9],
    Int64[7],
    Int64[8],
    Int64[9],
    Int64[],
    Int64[],
    Int64[]
]

wf.intOPs = [
    [601.0 1.0 602.0 0.0; 1202.0 1.0 1201.0 0.0],
    [801.0 1.0 802.0 0.0; 1201.0 1.0 1202.0 0.0; 1402.0 1.0 1401.0 0.0],
    [1001.0 1.0 1002.0 0.0; 1401.0 1.0 1402.0 0.0; 1602.0 1.0 1601.0 0.0],
    [1201.0 1.0 1202.0 0.0],
    [1401.0 1.0 1402.0 0.0],
    [1601.0 1.0 1602.0 0.0],
    zeros(0, 4),
    zeros(0, 4),
    zeros(0, 4)
]

# Test both versions: with and without explicit floris parameter

# Version 1: Explicit floris parameter (optimal)
ub1 = create_unified_buffers(wf, floris)
println("Created buffers with explicit floris parameter")

# Version 2: Default parameters (should work with improved default)
ub2 = create_unified_buffers(wf)
println("Created buffers with default parameters (50 rotor points)")

iT = 1
@assert !isempty(wf.intOPs[iT]) "wf.intOPs[$iT] is empty"

# Run once to warm up (test both buffer versions)
setUpTmpWFAndRun!(ub1, wf, set, floris, wind)
setUpTmpWFAndRun!(ub2, wf, set, floris, wind)

# Benchmark using optimal buffers (with explicit floris)
bench = @benchmark setUpTmpWFAndRun!(ub1, wf, set, floris, wind)

mean_time = mean(bench.times) / 1e6  # ms
allocs = mean(bench.memory) / 1024   # KiB
println("Benchmark setUpTmpWFAndRun!: $(round(mean_time, digits=3)) ms, $(round(allocs, digits=2)) KiB allocated")