# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause
 
using FLORIDyn, Test, Statistics, MAT, DataFrames

matlab_file   = "test/data/input_setUpTmpWFAndRun_2_steps.mat"
after_interpolateOPs_T_file = "test/data/after_interpolateOPs_T.mat"
vars = matread(matlab_file)
vars_after_interpolateOPs_T = matread(after_interpolateOPs_T_file)
wf_dict = vars["T"]
wf_dict_ref = vars_after_interpolateOPs_T["T"]

if !isdefined(Main, :TestHelpers)
    include("test_helpers.jl")
end
using .TestHelpers

before_interpolateOPs_T_file = "test/data/before_interpolateOPs_T.mat"
vars_before_interpolateOPs_T = matread(before_interpolateOPs_T_file)
wf_dict_03 = vars_before_interpolateOPs_T["T"]

after_interpolateOPs_T_file = "test/data/after_interpolateOPs_T.mat"
vars_after_interpolateOPs_T = matread(after_interpolateOPs_T_file)
wf_dict_02 = vars_after_interpolateOPs_T["T"]

@testset "interpolateOPs_basic" begin
    global wf, wf_ref_02
    wf = wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T
    # Create buffers for interpolateOPs!
    intOPs_buffers = [zeros(length(wf.dep[iT]), 4) for iT in 1:wf.nT]
    dist_buffer = zeros(wf.nOP)
    sorted_indices_buffer = zeros(Int, wf.nOP)
    wf.intOPs = interpolateOPs!(intOPs_buffers, wf, dist_buffer, sorted_indices_buffer)            # line 378ff in floridyn_cl.jl

    wf_ref_02 = wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T
    if !compare_windFarms(wf_ref_02, wf; detailed=false, tolerance=1e-6)
        @warn "WindFarm does not match reference after interpolateOPs"
        # compare_windFarms(wf_ref_02, wf; detailed=true, tolerance=1e-6)
        println("Calculated: ", wf.intOPs); println()
        println("Reference: ", wf_ref_02.intOPs)
    end
    @test compare_windFarms(wf_ref_02, wf; detailed=false, tolerance=1e-6)
end

@testset "interpolateOPs!_basic" begin
    global wf_no_alloc, wf_ref_02_copy
    
    # Start with the same initial state as the allocating version
    wf_no_alloc = wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T
    wf_ref_02_copy = wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T (reference)
    
    # Pre-allocate buffers for the non-allocating version
    intOPs = [zeros(length(wf_no_alloc.dep[iT]), 4) for iT in 1:wf_no_alloc.nT]
    dist_buffer = zeros(wf_no_alloc.nOP)
    sorted_indices_buffer = zeros(Int, wf_no_alloc.nOP)
    
    # Call the non-allocating version
    result = interpolateOPs!(intOPs, wf_no_alloc, dist_buffer, sorted_indices_buffer)
    
    # Assign the result to the wind farm object (same way as the allocating version)
    wf_no_alloc.intOPs = result
    
    # Test that the non-allocating version produces the same results as the reference
    if !compare_windFarms(wf_ref_02_copy, wf_no_alloc; detailed=false, tolerance=1e-6)
        @warn "WindFarm does not match reference after interpolateOPs!"
        # compare_windFarms(wf_ref_02_copy, wf_no_alloc; detailed=true, tolerance=1e-6)
        println("Calculated (non-allocating): ", wf_no_alloc.intOPs); println()
        println("Reference: ", wf_ref_02_copy.intOPs)
    end
    @test compare_windFarms(wf_ref_02_copy, wf_no_alloc; detailed=false, tolerance=1e-6)
end

@testset "interpolateOPs!_vs_interpolateOPs" begin
    # Test that both functions produce identical results
    global wf_alloc, wf_no_alloc_comp
    
    # Test with the same starting wind farm state
    wf_alloc = wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T
    wf_no_alloc_comp = wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T (copy)
    
    # Run allocating version (now using non-allocating version)
    intOPs_buffers_alloc = [zeros(length(wf_alloc.dep[iT]), 4) for iT in 1:wf_alloc.nT]
    dist_buffer_alloc = zeros(wf_alloc.nOP)
    sorted_indices_buffer_alloc = zeros(Int, wf_alloc.nOP)
    wf_alloc.intOPs = interpolateOPs!(intOPs_buffers_alloc, wf_alloc, dist_buffer_alloc, sorted_indices_buffer_alloc)
    
    # Run non-allocating version
    intOPs = [zeros(length(wf_no_alloc_comp.dep[iT]), 4) for iT in 1:wf_no_alloc_comp.nT]
    dist_buffer = zeros(wf_no_alloc_comp.nOP)
    sorted_indices_buffer = zeros(Int, wf_no_alloc_comp.nOP)
    
    wf_no_alloc_comp.intOPs = interpolateOPs!(intOPs, wf_no_alloc_comp, dist_buffer, sorted_indices_buffer)
    
    # Compare results between both versions
    @test compare_windFarms(wf_alloc, wf_no_alloc_comp; detailed=false, tolerance=1e-14)
    
    # Test that the specific intOPs field matches exactly
    @test length(wf_alloc.intOPs) == length(wf_no_alloc_comp.intOPs)
    for iT in 1:length(wf_alloc.intOPs)
        @test size(wf_alloc.intOPs[iT]) == size(wf_no_alloc_comp.intOPs[iT])
        @test wf_alloc.intOPs[iT] ≈ wf_no_alloc_comp.intOPs[iT] atol=1e-14
    end
end

@testset "interpolateOPs!_buffer_reuse" begin
    # Test that buffers can be safely reused across multiple calls
    global wf_buffer_test1, wf_buffer_test2
    
    # Create two different wind farm states
    wf_buffer_test1 = wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T
    wf_buffer_test2 = wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T (copy)
    
    # Pre-allocate buffers once
    max_nT = max(wf_buffer_test1.nT, wf_buffer_test2.nT)
    max_nOP = max(wf_buffer_test1.nOP, wf_buffer_test2.nOP)
    
    intOPs1 = [zeros(length(wf_buffer_test1.dep[iT]), 4) for iT in 1:wf_buffer_test1.nT]
    intOPs2 = [zeros(length(wf_buffer_test2.dep[iT]), 4) for iT in 1:wf_buffer_test2.nT]
    dist_buffer = zeros(max_nOP)
    sorted_indices_buffer = zeros(Int, max_nOP)
    
    # Run first calculation
    result1 = interpolateOPs!(intOPs1, wf_buffer_test1, dist_buffer, sorted_indices_buffer)
    wf_buffer_test1.intOPs = result1
    
    # Run second calculation with same buffers (simulating reuse)
    result2 = interpolateOPs!(intOPs2, wf_buffer_test2, dist_buffer, sorted_indices_buffer)
    wf_buffer_test2.intOPs = result2
    
    # Both should match the reference
    wf_ref_buffer1 = wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T
    wf_ref_buffer2 = wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T
    
    @test compare_windFarms(wf_ref_buffer1, wf_buffer_test1; detailed=false, tolerance=1e-6)
    @test compare_windFarms(wf_ref_buffer2, wf_buffer_test2; detailed=false, tolerance=1e-6)
    
    # Results should be identical since we started with the same state
    @test compare_windFarms(wf_buffer_test1, wf_buffer_test2; detailed=false, tolerance=1e-14)
end

nothing
