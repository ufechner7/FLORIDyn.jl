# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, BenchmarkTools

settings_file = "data/2021_9T_Data.yaml"
# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)
# create settings struct
set = Settings(wind, sim, con)
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf.dep = findTurbineGroups(wf, floridyn)
# Create unified buffers for interpolateOPs!
intOPs_buffers = [zeros(length(wf.dep[iT]), 4) for iT in 1:wf.nT]
unified_buffers = create_unified_buffers(wf)
wf.intOPs = interpolateOPs!(unified_buffers, intOPs_buffers, wf)
wf_old = deepcopy(wf)
@btime setUpTmpWFAndRun!(unified_buffers, wf, set, floris, wind)
nothing

# Laptop on battery
# 116.358 μs (2357 allocations: 397.69 KiB)
    