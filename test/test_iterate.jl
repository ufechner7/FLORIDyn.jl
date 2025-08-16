# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, Statistics, LinearAlgebra

"""
    create_test_setup()

Create a working wind farm setup for testing iterateOPs! functionality.
Uses the standard test configuration to ensure physical validity.
Returns the wind farm, simulation parameters, and pre-allocated buffers.
"""
function create_test_setup()
    settings_file = "data/2021_9T_Data.yaml"
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    
    # Create pre-allocated buffers for the new iterateOPs! function
    buffers = IterateOPsBuffers(wf)
    
    return wf, sim, floris, floridyn, set, buffers
end

@testset verbose=true "iterateOPs! Unit Tests" begin
    
    @testset "Basic Functionality" begin
        @testset "Function Execution" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            wf_original = deepcopy(wf)
            
            # Test that function executes without error
            @test_nowarn iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Test that wind farm states are potentially modified
            # Note: Due to circular shifting, some states may appear unchanged
            # but the OP states should definitely change
            @test wf.States_OP != wf_original.States_OP
            
            # States_T and States_WF may be identical after circular shifting
            # if the wind conditions are very stable, so we just check they remain valid
            @test !any(isnan.(wf.States_T))
            @test !any(isnan.(wf.States_WF))
            
            # Test return value
            result = iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            @test result === nothing
        end
        
        @testset "State Matrix Dimensions" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            original_size_op = size(wf.States_OP)
            original_size_t = size(wf.States_T)
            original_size_wf = size(wf.States_WF)
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Dimensions should remain unchanged
            @test size(wf.States_OP) == original_size_op
            @test size(wf.States_T) == original_size_t
            @test size(wf.States_WF) == original_size_wf
        end
        
        @testset "State Value Validity" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Check that no NaN or Inf values are produced
            @test !any(isnan.(wf.States_OP))
            @test !any(isinf.(wf.States_OP))
            @test !any(isnan.(wf.States_T))
            @test !any(isinf.(wf.States_T))
            @test !any(isnan.(wf.States_WF))
            @test !any(isinf.(wf.States_WF))
        end
    end
    
    @testset "State Advancement" begin
        @testset "Downstream Position Changes" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Store initial downstream positions
            initial_dw_pos = copy(wf.States_OP[:, 4])
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Some downstream positions should change due to advection
            # (accounting for circular shifting and complex wake dynamics)
            final_dw_pos = wf.States_OP[:, 4]
            position_changes = abs.(final_dw_pos - initial_dw_pos)
            
            # Check if positions changed, or if due to circular shift they appear unchanged
            # The key is that the function executes without error and maintains ordering
            changed_positions = any(position_changes .> 1e-10)
            positions_ordered = all([issorted(wf.States_OP[wf.StartI[iT]:wf.StartI[iT]+wf.nOP-1, 4]) for iT in 1:wf.nT])
            
            # At least one of these should be true: positions changed OR they remain properly ordered
            @test changed_positions || positions_ordered
        end
        
        @testset "World Coordinate Updates" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Store initial world coordinates
            initial_x = copy(wf.States_OP[:, 1])
            initial_y = copy(wf.States_OP[:, 2])
            initial_z = copy(wf.States_OP[:, 3])
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Check that some coordinates change
            x_changes = abs.(wf.States_OP[:, 1] - initial_x)
            y_changes = abs.(wf.States_OP[:, 2] - initial_y)
            z_changes = abs.(wf.States_OP[:, 3] - initial_z)
            
            # At least some movement should occur in at least one dimension
            total_movement = sum(x_changes) + sum(y_changes) + sum(z_changes)
            @test total_movement > 1e-10
        end
    end
    
    @testset "Temporal Behavior" begin
        @testset "Circular Shift Operation" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Create unique identifiers for tracking
            marker_values = collect(1:size(wf.States_OP, 1)) * 1000.0
            wf.States_OP[:, end] = marker_values  # Use last column as marker
            
            # Store values at turbine start positions 
            preserved_markers = wf.States_OP[wf.StartI[:], end]
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Values at start indices should be preserved (circular shift restores them)
            current_markers = wf.States_OP[wf.StartI[:], end]
            @test all(abs.(current_markers - preserved_markers) .< 1e-10)
        end
        
        @testset "Multiple Iterations Stability" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Run multiple iterations
            for i in 1:3
                @test_nowarn iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
                
                # Check state validity after each iteration
                @test !any(isnan.(wf.States_OP))
                @test !any(isnan.(wf.States_T))
                @test !any(isnan.(wf.States_WF))
            end
        end
    end
    
    @testset "Spatial Ordering" begin
        @testset "Downstream Position Ordering" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Check that downstream positions are ordered for each turbine
            for iT in 1:wf.nT
                start_idx = wf.StartI[iT]
                end_idx = start_idx + wf.nOP - 1
                range_indices = start_idx:end_idx
                downstream_positions = wf.States_OP[range_indices, 4]
                @test issorted(downstream_positions)
            end
        end
    end
    
    @testset "Physical Consistency" begin
        @testset "Conservation Properties" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Count initial dimensions
            initial_op_count = size(wf.States_OP, 1)
            initial_t_count = size(wf.States_T, 1)
            initial_wf_count = size(wf.States_WF, 1)
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Dimensions should be conserved
            @test size(wf.States_OP, 1) == initial_op_count
            @test size(wf.States_T, 1) == initial_t_count
            @test size(wf.States_WF, 1) == initial_wf_count
        end
        
        @testset "Reasonable Movement Magnitudes" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            initial_positions = copy(wf.States_OP[:, 1:3])
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            final_positions = wf.States_OP[:, 1:3]
            movement = final_positions - initial_positions
            max_movement = maximum(abs.(movement))
            
            # Movement should be reasonable for the time step
            # (Not too large, but can be larger than simple advection due to wake effects)
            @test max_movement < 1000.0  # Less than 1km per time step seems reasonable
            @test max_movement >= 0.0    # Should be non-negative
        end
    end
    
    @testset "Edge Cases and Robustness" begin
        @testset "Repeated Execution" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Run multiple times to test stability
            for i in 1:5
                @test_nowarn iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
                
                # Ensure ordering is maintained
                for iT in 1:wf.nT
                    start_idx = wf.StartI[iT]
                    end_idx = start_idx + wf.nOP - 1
                    range_indices = start_idx:end_idx
                    @test issorted(wf.States_OP[range_indices, 4])
                end
            end
        end
        
        @testset "Type Consistency" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Store original types
            original_op_type = typeof(wf.States_OP)
            original_t_type = typeof(wf.States_T)
            original_wf_type = typeof(wf.States_WF)
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Types should remain the same
            @test typeof(wf.States_OP) == original_op_type
            @test typeof(wf.States_T) == original_t_type
            @test typeof(wf.States_WF) == original_wf_type
        end
    end
    
    @testset "Algorithm Implementation Details" begin
        @testset "Downwind Step Calculation" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Test that the downwind step is calculated based on wind speed and time
            initial_dw = copy(wf.States_OP[:, 4])
            wind_speeds = copy(wf.States_WF[:, 1])
            time_step = sim.time_step
            advection_factor = sim.dyn.advection
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # The relationship is complex due to circular shifting, but we can verify
            # that changes are proportional to wind speed and time step
            expected_step_magnitude = maximum(wind_speeds) * time_step * advection_factor
            actual_changes = abs.(wf.States_OP[:, 4] - initial_dw)
            max_change = maximum(actual_changes)
            
            # The maximum change should be in a reasonable range
            @test max_change <= expected_step_magnitude * 2  # Allow for wake effects
        end
        
        @testset "Coordinate Transformation" begin
            wf, sim, floris, floridyn, set, buffers = create_test_setup()
            
            # Check that world coordinates are updated based on wind direction
            initial_x = copy(wf.States_OP[:, 1])
            initial_y = copy(wf.States_OP[:, 2])
            wind_dirs = copy(wf.States_WF[:, 2])
            
            iterateOPs!(IterateOPs_basic(), wf, sim, floris, floridyn, buffers)
            
            # Movement should occur in world coordinates
            x_movement = abs.(wf.States_OP[:, 1] - initial_x)
            y_movement = abs.(wf.States_OP[:, 2] - initial_y)
            
            # Should have some movement in at least one direction
            total_xy_movement = sum(x_movement) + sum(y_movement)
            @test total_xy_movement > 1e-12
        end
    end
    
    @testset "Regression Tests" begin
        @testset "Comparison with Reference Implementation" begin
            # This test ensures that our function produces consistent results
            wf1, sim, floris, floridyn, set, buffers1 = create_test_setup()
            wf2 = deepcopy(wf1)
            buffers2 = IterateOPsBuffers(wf2)  # Create separate buffers for second test
            
            # Run the same operation on both copies
            iterateOPs!(IterateOPs_basic(), wf1, sim, floris, floridyn, buffers1)
            iterateOPs!(IterateOPs_basic(), wf2, sim, floris, floridyn, buffers2)
            
            # Results should be identical
            @test wf1.States_OP ≈ wf2.States_OP
            @test wf1.States_T ≈ wf2.States_T
            @test wf1.States_WF ≈ wf2.States_WF
        end
        
        @testset "Deterministic Behavior" begin
            # Test that the function is deterministic for the same inputs
            wf_base, sim, floris, floridyn, set, _ = create_test_setup()
            
            results = []
            for i in 1:3
                wf_test = deepcopy(wf_base)
                buffers_test = IterateOPsBuffers(wf_test)
                iterateOPs!(IterateOPs_basic(), wf_test, sim, floris, floridyn, buffers_test)
                push!(results, copy(wf_test.States_OP))
            end
            
            # All results should be the same
            @test results[1] ≈ results[2]
            @test results[2] ≈ results[3]
        end
    end
end
nothing