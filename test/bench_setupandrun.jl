# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, BenchmarkTools, Parameters

import Base: show

settings_file = "data/2021_9T_Data.yaml"
wind, sim, con, floris, floridyn, ta = setup(settings_file)
set = Settings(wind, sim, con)

@with_kw mutable struct Allocs
    n::Int64 = 0      # number of floris calls
    m::Int64 = 0      # number of setupandrun calls
    for1::Int64 = 0   # allocated memory outer for loop
    for2::Int64 = 0   # allocated memory first inner for loop
    for3::Int64 = 0   # allocated memory first inner for loop
    if1::Int64 = 0    # allocated memory first if clause
    if2::Int64 = 0    # allocated memory first if clause
    if3::Int64 = 0    # allocated memory first if clause
    begin1::Int64 = 0
    floris::Int64 = 0
end

function Base.show(io::IO, allocs::Allocs)
    println(io, "Allocations:")
    for field_name in fieldnames(typeof(allocs))
        value = getfield(allocs, field_name)
        m = getfield(allocs, :m)
        if value > 5e2 && ! (field_name in [:n, :m])
            gb_value = value / 1024/ m
            println(io, "  $field_name: $(round(gb_value, digits=3)) KiB")
        end
    end
end

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

# Prepare unified buffers for setUpTmpWFAndRun!
ub = create_unified_buffers(wf)

alloc=Allocs()

iT = 1
@assert !isempty(wf.intOPs[iT]) "wf.intOPs[$iT] is empty"

# Run once to warm up
setUpTmpWFAndRun!(ub, wf, set, floris, wind)

# Benchmark
bench = @benchmark setUpTmpWFAndRun!(ub, wf, set, floris, wind; alloc)

mean_time = mean(bench.times) / 1e6  # ms
allocs = mean(bench.memory) / 1024   # KiB
println("Benchmark setUpTmpWFAndRun!: $(round(mean_time, digits=3)) ms, $(round(allocs, digits=2)) KiB allocated")
println(alloc)