# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, LinearAlgebra

"""
    create_test_setup_average()

Create a wind farm setup for testing iterateOPs! with IterateOPs_average.
Returns (wf, sim, floris, floridyn, set, buffers).
"""
function create_test_setup_average()
    settings_file = "data/2021_9T_Data.yaml"
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    buffers = IterateOPsBuffers(wf)
    return wf, sim, floris, floridyn, set, buffers
end

# Helper to set weights (mutate existing vector in-place)
function set_iter_weights!(sim::Sim, wNew::Float64, wOld::Float64)
    sim.dyn.op_iter_weights[1] = wNew
    sim.dyn.op_iter_weights[2] = wOld
    return sim
end

@testset "iterateOPs_average Unit Tests" begin

    @testset "Weight Identity (old state)" begin
        wf, sim, floris, floridyn, set, buffers = create_test_setup_average()
        # Store original wind field
        original_wf = copy(wf.States_WF)
        # Set weights so result should equal old (wNew=0, wOld=1)
        set_iter_weights!(sim, 0.0, 1.0)
        @test_nowarn iterateOPs!(IterateOPs_average(), wf, sim, floris, floridyn, buffers)
        @test wf.States_WF ≈ original_wf atol=1e-12
    end

    @testset "Weight Identity (shifted state matches basic)" begin
        # Two identical setups
        wf_a, sim_a, floris_a, floridyn_a, set_a, buffers_a = create_test_setup_average()
        wf_b = deepcopy(wf_a); sim_b = deepcopy(sim_a); floris_b = deepcopy(floris_a); floridyn_b = floridyn_a; buffers_b = IterateOPsBuffers(wf_b)
        # Configure weights wNew=1 (shifted), wOld=0
        set_iter_weights!(sim_a, 1.0, 0.0)
        iterateOPs!(IterateOPs_average(), wf_a, sim_a, floris_a, floridyn_a, buffers_a)
        iterateOPs!(IterateOPs_basic(),   wf_b, sim_b, floris_b, floridyn_b, buffers_b)
        @test wf_a.States_WF ≈ wf_b.States_WF atol=1e-12
    end

    @testset "General Weighted Combination" begin
        wf, sim, floris, floridyn, set, buffers = create_test_setup_average()
        # Capture original state
        old = copy(wf.States_WF)
        # Build initial_states matrix for start rows (as done in iterate)
        nT = size(wf.StartI,2)
        initial_states = similar(buffers.tmpWFStates)
        for i in 1:nT
            start_idx = wf.StartI[1,i]
            @views initial_states[i,:] .= old[start_idx, :]
        end
        # Prepare a copy to produce the shifted+restored matrix S
        shifted = copy(old)
        temp_buffer = similar(shifted)
        # Use internal helper (not exported)
        FLORIDyn._circshift_and_restore!(shifted, initial_states, wf.StartI, temp_buffer)
        # Choose non-trivial weights
        wNew = 0.3; wOld = 0.7
        set_iter_weights!(sim, wNew, wOld)
        # Expected following Julia implementation order:
        expected = wOld .* old .+ wNew .* shifted
        # Restore start rows (Julia restores after averaging)
        for i in 1:nT
            start_idx = wf.StartI[1,i]
            @views expected[start_idx, :] .= old[start_idx, :]
        end
        iterateOPs!(IterateOPs_average(), wf, sim, floris, floridyn, buffers)
        @test wf.States_WF ≈ expected rtol=1e-12 atol=1e-12
    end

    @testset "Deterministic Behavior" begin
        wf, sim, floris, floridyn, set, buffers = create_test_setup_average()
        set_iter_weights!(sim, 0.4, 0.6)
        results = Matrix{Float64}[]
        for i in 1:3
            wf_copy = deepcopy(wf)
            sim_copy = deepcopy(sim)
            buffers_copy = IterateOPsBuffers(wf_copy)
            iterateOPs!(IterateOPs_average(), wf_copy, sim_copy, floris, floridyn, buffers_copy)
            push!(results, copy(wf_copy.States_WF))
        end
        @test results[1] ≈ results[2]
        @test results[2] ≈ results[3]
    end
end
nothing
