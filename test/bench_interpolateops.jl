# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, Statistics, LinearAlgebra, BenchmarkTools

"""
    create_test_setup()

Create a working wind farm setup for testing iterateOPs! functionality.
Uses the standard test configuration to ensure physical validity.
"""
function create_test_setup()
    settings_file = "data/2021_9T_Data.yaml"
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    set = Settings(wind, sim, con)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    
    return wf, wind, sim, con, floris, floridyn, ta, set
end


wf, wind, sim, con, floris, floridyn, ta, set = create_test_setup()
wf_original = deepcopy(wf)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf.dep = findTurbineGroups(wf, floridyn)
# Create unified buffers for interpolateOPs!
intOPs_buffers = [zeros(length(wf.dep[iT]), 4) for iT in 1:wf.nT]
unified_buffers = create_unified_buffers(wf)
t = @benchmark wf.intOPs = interpolateOPs!($intOPs_buffers, $wf, $unified_buffers)

time = mean(t.times)/1e9
rel_time = time * 301 / 0.08  # Relative to the total time of 0.08 seconds
println("Benchmark time: $time seconds, relative to 0.08s: $(round(rel_time * 100, digits=2)) %")

# Benchmark time: 1.99795341e-5 seconds, relative to 0.115s: 5.23 %