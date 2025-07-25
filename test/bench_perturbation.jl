# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, BenchmarkTools

settings_file = "data/2021_9T_Data.yaml"
# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn = setup(settings_file)
# create settings struct
set = Settings(wind, sim, con)
# % Load linked data
turbProp        = turbineArrayProperties(settings_file)
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbProp, sim)
wf_old = deepcopy(wf)
@btime perturbationOfTheWF!(wf, wind)
wind.pertubation.vel = true
@btime perturbationOfTheWF!(wf, wind)
wind.pertubation.vel = false
wind.pertubation.dir = true
perturbationOfTheWF!(wf, wind)
wind.pertubation.dir = false
wind.pertubation.ti = true
@btime perturbationOfTheWF!(wf, wind)

# On AMD 7890X:
# 9.830 ns (0 allocations: 0 bytes)
# 2.681 μs (9 allocations: 42.40 KiB)
# 2.667 μs (9 allocations: 42.40 KiB)