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
@btime vv_dep = findTurbineGroups(wf, floridyn)
nothing

# On AMD 7890X:
# 97.180 μs (1874 allocations: 626.83 KiB)