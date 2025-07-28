# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, BenchmarkTools

settings_file = "data/2021_9T_Data.yaml"
# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)
# create settings struct
set = Settings(wind, sim, con)
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf_old = deepcopy(wf)
buffers = FLORIDyn.IterateOPsBuffers(wf)
t = @benchmark iterateOPs!(set.iterate_mode, wf, sim, floris, floridyn, buffers)

time = mean(t.times)/1e9
rel_time = time * 301 / 0.115  # Relative to the total time of 0.115 seconds
println("Benchmark time: $time seconds, relative to 0.115s: $(round(rel_time * 100, digits=2)) %")

# On AMD 7890X:
# 149 µs (180 allocations: 892.77 KiB)
# 136 μs  (90 allocations: 579.94 KiB) (using views)
#  82 μs  (0 allocations)  (allocation free version with preallocated buffers)
# Benchmark time: 0.00010247352280000001 seconds, relative to 0.115s: 26.82 %
# Benchmark time: 8.33199938e-5 seconds, relative to 0.115s: 21.81 %