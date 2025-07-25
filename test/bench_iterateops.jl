# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, BenchmarkTools

settings_file = "data/2021_9T_Data.yaml"
# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn = setup(settings_file)
# create settings struct
set = Settings(wind, sim, con)
# % Load linked data
turbine_prop        = turbineArrayProperties(settings_file)
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_prop, sim)
wf_old = deepcopy(wf)
@btime iterateOPs!(set.iterate_mode, wf, sim, floris, floridyn)

# 149 µs (on AMD 7890X)
