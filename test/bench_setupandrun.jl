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
        if value > 5e2 && field_name != :m
            gb_value = value / 1024/ m
            println(io, "  $field_name: $(round(gb_value, digits=3)) KiB")
        end
    end
end

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf = initSimulation(wf, sim)

# Initialize dependencies for all turbines
wf.dep = [Int64[] for _ in 1:wf.nT]

# Prepare buffers for setUpTmpWFAndRun!
nT = wf.nT
nOP = wf.nOP
nWF = size(wf.States_WF, 2)
nTst = size(wf.States_T, 2)

M_buffer = zeros(nT, 3)
iTWFState_buffer = zeros(nWF)
tmp_Tpos_buffer = zeros(nT, 3)
tmp_WF_buffer = zeros(nT, nWF)
tmp_Tst_buffer = zeros(nT, nTst)
dists_buffer = zeros(nT)
plot_WF_buffer = zeros(nT, nWF)
plot_OP_buffer = zeros(nT, 2)

alloc=Allocs()

# Run once to warm up
setUpTmpWFAndRun!(M_buffer, wf, set, floris, wind,
    iTWFState_buffer, tmp_Tpos_buffer, tmp_WF_buffer, tmp_Tst_buffer,
    dists_buffer, plot_WF_buffer, plot_OP_buffer)

# Benchmark
bench = @benchmark setUpTmpWFAndRun!(M_buffer, wf, set, floris, wind,
    iTWFState_buffer, tmp_Tpos_buffer, tmp_WF_buffer, tmp_Tst_buffer,
    dists_buffer, plot_WF_buffer, plot_OP_buffer; alloc)

mean_time = mean(bench.times) / 1e6  # ms
allocs = mean(bench.memory) / 1024   # KiB
println("Benchmark setUpTmpWFAndRun!: $(round(mean_time, digits=3)) ms, $(round(allocs, digits=2)) KiB allocated")
println(alloc)