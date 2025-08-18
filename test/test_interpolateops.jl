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

@testset verbose = true "InterpolateOPs Tests" begin
    @testset "interpolateOPs_basic" begin
        global wf, wf_ref_02
        wf = wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T
        # Create unified buffers for interpolateOPs!
        unified_buffers = create_unified_buffers(wf)
        wf.intOPs = [zeros(length(wf.dep[iT]), 4) for iT in 1:wf.nT]
        interpolateOPs!(unified_buffers, wf.intOPs, wf)

        wf_ref_02 = wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T
        if !compare_windFarms(wf_ref_02, wf; detailed=false, tolerance=1e-6)
            @warn "WindFarm does not match reference after interpolateOPs"
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
        result = [zeros(length(wf_no_alloc.dep[iT]), 4) for iT in 1:wf_no_alloc.nT]
        unified_buffers = create_unified_buffers(wf_no_alloc)
        
        # Call the non-allocating version
        interpolateOPs!(unified_buffers, result, wf_no_alloc)

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
        wf_alloc.intOPs = [zeros(length(wf_alloc.dep[iT]), 4) for iT in 1:wf_alloc.nT]
        unified_buffers_alloc = create_unified_buffers(wf_alloc)
        interpolateOPs!(unified_buffers_alloc, wf_alloc.intOPs, wf_alloc)

        # Run non-allocating version
        wf_no_alloc_comp.intOPs = [zeros(length(wf_no_alloc_comp.dep[iT]), 4) for iT in 1:wf_no_alloc_comp.nT]
        unified_buffers = create_unified_buffers(wf_no_alloc_comp)

        interpolateOPs!(unified_buffers, wf_no_alloc_comp.intOPs, wf_no_alloc_comp)

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
        
        result1 = [zeros(length(wf_buffer_test1.dep[iT]), 4) for iT in 1:wf_buffer_test1.nT]
        result2 = [zeros(length(wf_buffer_test2.dep[iT]), 4) for iT in 1:wf_buffer_test2.nT]
        unified_buffers = create_unified_buffers(wf_buffer_test1)  # Both use same wind farm structure
        
        # Run first calculation
        interpolateOPs!(unified_buffers, result1, wf_buffer_test1)
        wf_buffer_test1.intOPs = result1
        
        # Run second calculation with same buffers (simulating reuse)
        interpolateOPs!(unified_buffers, result2, wf_buffer_test2)
        wf_buffer_test2.intOPs = result2
        
        # Both should match the reference
        wf_ref_buffer1 = wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T
        wf_ref_buffer2 = wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T
        
        @test compare_windFarms(wf_ref_buffer1, wf_buffer_test1; detailed=false, tolerance=1e-6)
        @test compare_windFarms(wf_ref_buffer2, wf_buffer_test2; detailed=false, tolerance=1e-6)
        
        # Results should be identical since we started with the same state
        @test compare_windFarms(wf_buffer_test1, wf_buffer_test2; detailed=false, tolerance=1e-14)
    end
end
nothing
