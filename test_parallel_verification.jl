# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Test script to verify parallel execution is working correctly
using FLORIDyn
using ControlPlots

function test_parallel_vs_sequential()
    println("Testing parallel vs sequential execution...")
    
    # Setup common parameters
    wind, sim, con, floris, floridyn, ta = setup("data/2021_9T_Data.yaml")
    
    # Test sequential execution
    println("Testing sequential execution...")
    set_seq = Settings(wind, sim, con)
    set_seq.parallel = false
    wf_seq, wind_seq, sim_seq, con_seq, floris_seq = prepareSimulation(set_seq, wind, con, floridyn, floris, ta, sim)
    wf_seq = initSimulation(wf_seq, sim_seq)
    
    t_seq = @elapsed Z_seq, X_seq, Y_seq = calcFlowField(set_seq, wf_seq, wind_seq, floris_seq)
    println("Sequential time: $(t_seq) seconds")
    
    # Test parallel execution  
    println("Testing parallel execution...")
    set_par = Settings(wind, sim, con)
    set_par.parallel = true
    wf_par, wind_par, sim_par, con_par, floris_par = prepareSimulation(set_par, wind, con, floridyn, floris, ta, sim)
    wf_par = initSimulation(wf_par, sim_par)
    
    t_par = @elapsed Z_par, X_par, Y_par = calcFlowField(set_par, wf_par, wind_par, floris_par)
    println("Parallel time: $(t_par) seconds")
    
    # Verify results are similar (they should be identical or very close)
    println("Checking results consistency...")
    max_diff = maximum(abs.(Z_seq - Z_par))
    println("Maximum difference between sequential and parallel results: $max_diff")
    
    if max_diff < 1e-10
        println("✓ Results are identical!")
    elseif max_diff < 1e-6
        println("✓ Results are very close (within numerical precision)")
    else
        println("⚠ Results differ significantly!")
        return false
    end
    
    speedup = t_seq / t_par
    println("Speedup: $(round(speedup, digits=2))x")
    
    return true
end

# Run the test
test_parallel_vs_sequential()
