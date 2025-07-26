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
    wind, sim, con, floris, floridyn = setup(settings_file)
    set = Settings(wind, sim, con)
    turbProp = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbProp, sim)
    
    return wf, sim, floris, floridyn, set
end

    
wf, sim, floris, floridyn, set = create_test_setup()
wf_original = deepcopy(wf)

# Test that function executes without error
t = @benchmark iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn)

time = mean(t.times)/1e9
rel_time = time * 301 / 0.115  # Relative to the total time of 0.115 seconds
println("Benchmark time: $time seconds, relative to 0.115s: $(round(rel_time * 100, digits=2)) %")

# Benchmark time: 0.0001022896801 seconds, relative to 0.115s: 26.77 %