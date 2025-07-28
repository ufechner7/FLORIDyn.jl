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
@btime perturbationOfTheWF!(wf, wind)
wind.perturbation.vel = true
@btime perturbationOfTheWF!(wf, wind)
wind.perturbation.vel = false
wind.perturbation.dir = true
perturbationOfTheWF!(wf, wind)
wind.perturbation.dir = false
wind.perturbation.ti = true
t = @benchmark perturbationOfTheWF!(wf, wind)

time = mean(t.times)/1e9
rel_time = time * 301 / 0.115  # Relative to the total time of 0.115 seconds
println("Benchmark time: $time seconds, relative to 0.115s: $(round(rel_time * 100, digits=2)) %")

# On AMD 7890X:
# 9.830 ns (0 allocations: 0 bytes)
# 2.681 μs (9 allocations: 42.40 KiB)
# 2.667 μs (9 allocations: 42.40 KiB)
# Benchmark time: 3.1538405666666664e-6 seconds, relative to 0.115s: 0.83 %