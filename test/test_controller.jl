# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, Interpolations

# Import the getYaw function directly
import FLORIDyn: getYaw

@testset verbose=true "getYaw function tests" begin
    
    @testset "Basic interpolation with single turbine" begin
        # Test data: time column + one turbine yaw column
        ConYawData = [
            0.0  10.0;
            1.0  15.0;
            2.0  20.0;
            3.0  25.0;
            4.0  30.0
        ]
        con=Con(yaw="Interpolate")
        con.yaw_data=ConYawData
        
        # Test interpolation at exact time points
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 0.0) ≈ 10.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 1.0) ≈ 15.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 2.0) ≈ 20.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 3.0) ≈ 25.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 4.0) ≈ 30.0
        
        # Test interpolation at intermediate time points
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 0.5) ≈ 12.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 1.5) ≈ 17.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 2.5) ≈ 22.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 3.5) ≈ 27.5
    end
    
    @testset "Multiple turbines interpolation" begin
        # Test data: time column + three turbine yaw columns
        ConYawData = [
            0.0  10.0  5.0   0.0;
            1.0  15.0  10.0  5.0;
            2.0  20.0  15.0  10.0;
            3.0  25.0  20.0  15.0
        ]

        con=Con(yaw="Interpolate")
        con.yaw_data=ConYawData
        
        # Test single turbine access
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 1.5) ≈ 17.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 2, 1.5) ≈ 12.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 3, 1.5) ≈ 7.5
        
        # Test multiple turbine access
        result = getYaw(FLORIDyn.Yaw_SOWFA(), con, [1, 2, 3], 1.5)
        @test result ≈ [17.5, 12.5, 7.5]
        
        # Test subset of turbines
        result = getYaw(FLORIDyn.Yaw_SOWFA(), con, [1, 3], 2.5)
        @test result ≈ [22.5, 12.5]
        
        # Test single turbine in vector format
        result = getYaw(FLORIDyn.Yaw_SOWFA(), con, [2], 0.5)
        @test result ≈ [7.5]
    end
    
    @testset "Out of bounds time handling" begin
        ConYawData = [
            1.0  10.0;
            2.0  20.0;
            3.0  30.0;
            4.0  40.0
        ]

        con = Con(yaw="SOWFA")
        con.yaw_data=ConYawData
        
        # Test time before first time point (should clamp to first time and issue warning)
        @test_logs (:warn, "The time 0.5 is out of bounds, will use 1.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 0.5)
            @test result ≈ 10.0
        end
        
        # Test time after last time point (should clamp to last time and issue warning)
        @test_logs (:warn, "The time 5.0 is out of bounds, will use 4.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 5.0)
            @test result ≈ 40.0
        end
        
        # Test with multiple turbines and out of bounds time
        ConYawData_multi = [
            1.0  10.0  5.0;
            2.0  20.0  15.0;
            3.0  30.0  25.0
        ]
        con.yaw_data=ConYawData_multi
        
        @test_logs (:warn, "The time 0.0 is out of bounds, will use 1.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), con, [1, 2], 0.0)
            @test result ≈ [10.0, 5.0]
        end
    end
    
    @testset "Edge cases and error handling" begin
        ConYawData = [
            0.0  10.0  5.0;
            1.0  20.0  15.0;
            2.0  30.0  25.0
        ]

        con = Con(yaw="SOWFA")
        con.yaw_data=ConYawData
        
        # Test invalid turbine index types
        @test_throws ErrorException getYaw(FLORIDyn.Yaw_SOWFA(), con, 1.5, 1.0)
        @test_throws ErrorException getYaw(FLORIDyn.Yaw_SOWFA(), con, "1", 1.0)
        @test_throws ErrorException getYaw(FLORIDyn.Yaw_SOWFA(), con, [1.5, 2.0], 1.0)
        
        # Test valid edge cases
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, Int32(1), 1.0) ≈ 20.0  # Different integer type
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, [Int32(1), Int64(2)], 1.0) ≈ [20.0, 15.0]  # Mixed integer types
    end
    
    @testset "Single time point data" begin
        # Test with only one time point
        ConYawData = [5.0  25.0  35.0]

        con = Con(yaw="SOWFA")
        con.yaw_data=ConYawData
        
        # Should return the single value regardless of requested time
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 5.0) ≈ 25.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 2, 5.0) ≈ 35.0
        
        # Out of bounds should still issue warnings but return the single value
        @test_logs (:warn, "The time 0.0 is out of bounds, will use 5.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 0.0)
            @test result ≈ 25.0
        end
        
        @test_logs (:warn, "The time 10.0 is out of bounds, will use 5.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), con, [1, 2], 10.0)
            @test result ≈ [25.0, 35.0]
        end
    end
    
    @testset "Two time point data" begin
        # Test with exactly two time points (minimum for linear interpolation)
        ConYawData = [
            0.0  0.0   10.0;
            10.0 100.0 20.0
        ]

        con = Con(yaw="SOWFA")
        con.yaw_data=ConYawData
        
        # Test exact endpoints
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 0.0) ≈ 0.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 10.0) ≈ 100.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 2, 0.0) ≈ 10.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 2, 10.0) ≈ 20.0
        
        # Test interpolation
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 5.0) ≈ 50.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 2, 5.0) ≈ 15.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, [1, 2], 2.5) ≈ [25.0, 12.5]
    end
    
    @testset "Realistic wind farm scenario" begin
        # Simulate a realistic scenario with multiple turbines and time-varying yaw
        # Time from 0 to 600 seconds (10 minutes), 9 turbines
        time_points = [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0]
        ConYawData = [
            0.0   0.0   0.0   0.0   5.0   5.0   5.0   10.0  10.0  10.0;
            100.0 0.0   5.0   10.0  5.0   10.0  15.0  10.0  15.0  20.0;
            200.0 5.0   10.0  15.0  10.0  15.0  20.0  15.0  20.0  25.0;
            300.0 10.0  15.0  20.0  15.0  20.0  25.0  20.0  25.0  30.0;
            400.0 5.0   10.0  15.0  10.0  15.0  20.0  15.0  20.0  25.0;
            500.0 0.0   5.0   10.0  5.0   10.0  15.0  10.0  15.0  20.0;
            600.0 0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0   0.0
        ]

        con = Con(yaw="SOWFA")
        con.yaw_data=ConYawData
        
        # Test at specific time points
        result_t150 = getYaw(FLORIDyn.Yaw_SOWFA(), con, 1:9, 150.0)
        expected_t150 = [2.5, 7.5, 12.5, 7.5, 12.5, 17.5, 12.5, 17.5, 22.5]
        @test result_t150 ≈ expected_t150
        
        # Test subset of turbines at different time
        result_subset = getYaw(FLORIDyn.Yaw_SOWFA(), con, [1, 5, 9], 250.0)
        expected_subset = [7.5, 17.5, 27.5]
        @test result_subset ≈ expected_subset
        
        # Test single turbine across different times
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 5, 0.0) ≈ 5.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 5, 300.0) ≈ 20.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), con, 5, 600.0) ≈ 0.0
    end
    
    @testset "Extrapolation behavior (Flat boundary condition)" begin
        # The function uses Flat() extrapolation, so values beyond bounds should be constant
        ConYawData = [
            1.0  10.0;
            2.0  20.0;
            3.0  30.0
        ]

        con = Con(yaw="SOWFA")
        con.yaw_data=ConYawData
        
        # Note: The function clamps time to bounds and issues warnings,
        # but the underlying interpolation would use flat extrapolation
        # We test this by checking the clamping behavior
        
        @test_logs (:warn, "The time 0.5 is out of bounds, will use 1.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 0.5)
            @test result ≈ 10.0  # Should be the first value
        end
        
        @test_logs (:warn, "The time 4.0 is out of bounds, will use 3.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), con, 1, 4.0)
            @test result ≈ 30.0  # Should be the last value
        end
    end
    
    @testset "Performance with large datasets" begin
        # Test with a larger dataset to ensure reasonable performance
        n_times = 1000
        n_turbines = 50
        
        times = collect(0.0:1.0:(n_times-1))
        yaw_angles = rand(n_times, n_turbines) * 360  # Random yaw angles 0-360 degrees
        
        ConYawData = hcat(times, yaw_angles)
        con = Con(yaw="SOWFA")
        con.yaw_data=ConYawData
        
        # Test single turbine access
        @time result1 = getYaw(FLORIDyn.Yaw_SOWFA(), con, 25, 500.5)
        @test isa(result1, Float64)
        
        # Test multiple turbine access
        @time result_all = getYaw(FLORIDyn.Yaw_SOWFA(), con, 1:n_turbines, 500.5)
        @test length(result_all) == n_turbines
        @test all(isa.(result_all, Float64))
        
        # Test subset access
        @time result_subset = getYaw(FLORIDyn.Yaw_SOWFA(), con, [1, 10, 25, 40, 50], 750.3)
        @test length(result_subset) == 5
    end
    
    @testset "Yaw_Constant tests" begin
        @testset "Basic functionality" begin
            # Test with single constant value
            con = Con(yaw="Constant")
            con.yaw_fixed = 225.0
            
            # Test single turbine
            @test getYaw(FLORIDyn.Yaw_Constant(), con, 1, 0.0) ≈ 225.0
            @test getYaw(FLORIDyn.Yaw_Constant(), con, 1, 100.0) ≈ 225.0  # Time ignored
            @test getYaw(FLORIDyn.Yaw_Constant(), con, 5, 500.0) ≈ 225.0  # Turbine index ignored
            
            # Test multiple turbines
            result = getYaw(FLORIDyn.Yaw_Constant(), con, [1, 2, 3], 0.0)
            @test result ≈ [225.0, 225.0, 225.0]
            
            # Test with larger vector
            result_large = getYaw(FLORIDyn.Yaw_Constant(), con, 1:10, 42.5)
            @test result_large ≈ fill(225.0, 10)
        end
        
        @testset "Error handling" begin
            # Test invalid turbine index types
            con = Con(yaw="Constant")
            con.yaw_fixed = 135.0
            @test_throws ErrorException getYaw(FLORIDyn.Yaw_Constant(), con, 1.5, 0.0)
            @test_throws ErrorException getYaw(FLORIDyn.Yaw_Constant(), con, "1", 0.0)
            @test_throws ErrorException getYaw(FLORIDyn.Yaw_Constant(), con, [1.5, 2.0], 0.0)
        end
        
        @testset "Edge cases and type variations" begin
            con = Con(yaw="Constant")
            con.yaw_fixed = 270.5
            
            # Test different integer types for turbine indices
            @test getYaw(FLORIDyn.Yaw_Constant(), con, Int32(1), 0.0) ≈ 270.5
            @test getYaw(FLORIDyn.Yaw_Constant(), con, Int64(1), 0.0) ≈ 270.5
            @test getYaw(FLORIDyn.Yaw_Constant(), con, UInt8(1), 0.0) ≈ 270.5
            
            # Test mixed integer types in vector
            result = getYaw(FLORIDyn.Yaw_Constant(), con, [Int32(1), Int64(2), UInt8(3)], 0.0)
            @test result ≈ [270.5, 270.5, 270.5]
            
            # Test with single element vector
            result_single = getYaw(FLORIDyn.Yaw_Constant(), con, [1], 0.0)
            @test result_single ≈ [270.5]
            @test isa(result_single, Vector{Float64})
        end
        
#         @testset "Real-world scenarios" begin
#             # Test with typical wind farm yaw angles
#             ConYawData = [0.0;;]  # Aligned with wind
#             result_aligned = getYaw(FLORIDyn.Yaw_Constant(), ConYawData, 1:54, 0.0)  # 54 turbines
#             @test length(result_aligned) == 54
#             @test all(result_aligned .≈ 0.0)
            
#             # Test with wake steering angle
#             ConYawData = [25.0;;]  # 25 degree wake steering
#             result_steering = getYaw(FLORIDyn.Yaw_Constant(), ConYawData, [1, 5, 10, 20], 1000.0)
#             @test result_steering ≈ [25.0, 25.0, 25.0, 25.0]
            
#             # Test with negative yaw angle
#             ConYawData = [-15.0;;]
#             result_negative = getYaw(FLORIDyn.Yaw_Constant(), ConYawData, 1, 0.0)
#             @test result_negative ≈ -15.0
            
#             # Test with large yaw angle
#             ConYawData = [359.9;;]
#             result_large_angle = getYaw(FLORIDyn.Yaw_Constant(), ConYawData, [1, 2, 3], 0.0)
#             @test result_large_angle ≈ [359.9, 359.9, 359.9]
#         end
        
#         @testset "Performance and consistency" begin
#             # Test that time parameter is truly ignored
#             ConYawData = [180.0;;]
            
#             times_to_test = [-1000.0, -1.0, 0.0, 1.0, 100.0, 1e6]
#             for t in times_to_test
#                 @test getYaw(FLORIDyn.Yaw_Constant(), ConYawData, 1, t) ≈ 180.0
#                 @test getYaw(FLORIDyn.Yaw_Constant(), ConYawData, [1, 2, 3], t) ≈ [180.0, 180.0, 180.0]
#             end
            
#             # Test with large number of turbines for performance
#             large_turbine_indices = 1:1000
#             @time result_perf = getYaw(FLORIDyn.Yaw_Constant(), ConYawData, large_turbine_indices, 0.0)
#             @test length(result_perf) == 1000
#             @test all(result_perf .≈ 180.0)
#         end
    end
# end

# # Import the getInduction function directly
# import FLORIDyn: getInduction

# @testset verbose=true "getInduction function tests" begin
    
#     @testset "Induction_Constant tests" begin
#         @testset "Basic constant induction functionality" begin
#             # Test data with single constant value
#             con_induction_data = [0.33;;]
            
#             # Test single turbine
#             @test getInduction(FLORIDyn.Induction_Constant(), con_induction_data, 1, 0.0) ≈ 0.33
#             @test getInduction(FLORIDyn.Induction_Constant(), con_induction_data, 1, 100.0) ≈ 0.33
#             @test getInduction(FLORIDyn.Induction_Constant(), con_induction_data, 1, 1000.0) ≈ 0.33
            
#             # Test multiple turbines
#             result = getInduction(FLORIDyn.Induction_Constant(), con_induction_data, [1, 2, 3], 50.0)
#             @test result == [0.33, 0.33, 0.33]
            
#             # Test with different matrix sizes (should ignore extra elements)
#             con_induction_data_large = [0.33 0.25 0.40; 
#                                        0.35 0.30 0.45]
#             @test getInduction(FLORIDyn.Induction_Constant(), con_induction_data_large, 1, 25.0) ≈ 0.33
#             @test getInduction(FLORIDyn.Induction_Constant(), con_induction_data_large, [1, 2], 75.0) == [0.33, 0.33]
#         end
        
#         @testset "Error handling for Induction_Constant" begin
#             # Test empty matrix error handling
#             empty_data_rows = Array{Float64}(undef, 0, 1)
#             @test_throws ErrorException getInduction(FLORIDyn.Induction_Constant(), empty_data_rows, 1, 0.0)
            
#             empty_data_cols = Array{Float64}(undef, 1, 0)
#             @test_throws ErrorException getInduction(FLORIDyn.Induction_Constant(), empty_data_cols, 1, 0.0)
            
#             # Test invalid iT parameter types
#             con_induction_data = [0.33;;]
#             @test_throws ErrorException getInduction(FLORIDyn.Induction_Constant(), con_induction_data, "invalid", 0.0)
#             @test_throws ErrorException getInduction(FLORIDyn.Induction_Constant(), con_induction_data, 1.5, 0.0)
#             @test_throws ErrorException getInduction(FLORIDyn.Induction_Constant(), con_induction_data, [1.5, 2.5], 0.0)
#         end
#     end
    
#     @testset "Induction_MPC tests" begin
#         @testset "Basic interpolation with single turbine" begin
#             # Test data: time column + one turbine induction column
#             con_induction_data = [
#                 0.0  0.30;
#                 1.0  0.32;
#                 2.0  0.34;
#                 3.0  0.36;
#                 4.0  0.38
#             ]
            
#             # Test interpolation at exact time points
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 0.0) ≈ 0.30
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 1.0) ≈ 0.32
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 2.0) ≈ 0.34
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 3.0) ≈ 0.36
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 4.0) ≈ 0.38
            
#             # Test interpolation at intermediate time points
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 0.5) ≈ 0.31
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 1.5) ≈ 0.33
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 2.5) ≈ 0.35
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 3.5) ≈ 0.37
#         end
        
#         @testset "Multiple turbines interpolation" begin
#             # Test data: time column + three turbine induction columns
#             con_induction_data = [
#                 0.0  0.30  0.25  0.20;
#                 1.0  0.32  0.27  0.22;
#                 2.0  0.34  0.29  0.24;
#                 3.0  0.36  0.31  0.26
#             ]
            
#             # Test single turbine access
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 1.5) ≈ 0.33
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 2, 1.5) ≈ 0.28
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 3, 1.5) ≈ 0.23
            
#             # Test multiple turbine access
#             result = getInduction(FLORIDyn.Induction_MPC(), con_induction_data, [1, 2, 3], 1.5)
#             @test result ≈ [0.33, 0.28, 0.23]
            
#             # Test partial turbine selection
#             result = getInduction(FLORIDyn.Induction_MPC(), con_induction_data, [2, 3], 2.5)
#             @test result ≈ [0.30, 0.25]
#         end
        
#         @testset "Out of bounds time handling" begin
#             con_induction_data = [
#                 0.0  0.30;
#                 1.0  0.32;
#                 2.0  0.34
#             ]
            
#             # Test time before range (should warn and clamp to first time)
#             @test_logs (:warn, r"The time -1.0 is out of bounds") begin
#                 result = getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, -1.0)
#                 @test result ≈ 0.30
#             end
            
#             # Test time after range (should warn and clamp to last time)
#             @test_logs (:warn, r"The time 5.0 is out of bounds") begin
#                 result = getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 5.0)
#                 @test result ≈ 0.34
#             end
#         end
        
#         @testset "Single time point data" begin
#             # Test data with only one time point
#             con_induction_data = [1.0  0.33  0.28  0.31]
            
#             # Test single turbine
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 0.0) ≈ 0.33
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 2, 5.0) ≈ 0.28
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 3, 1.0) ≈ 0.31
            
#             # Test multiple turbines
#             result = getInduction(FLORIDyn.Induction_MPC(), con_induction_data, [1, 2, 3], 2.0)
#             @test result ≈ [0.33, 0.28, 0.31]
#         end
        
#         @testset "Edge cases and error handling" begin
#             con_induction_data = [
#                 0.0  0.30  0.25;
#                 1.0  0.32  0.27;
#                 2.0  0.34  0.29
#             ]
            
#             # Test invalid iT parameter types
#             @test_throws ErrorException getInduction(FLORIDyn.Induction_MPC(), con_induction_data, "invalid", 1.0)
#             @test_throws ErrorException getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1.5, 1.0)
#             @test_throws ErrorException getInduction(FLORIDyn.Induction_MPC(), con_induction_data, [1.5, 2.5], 1.0)
#         end
        
#         @testset "Realistic wind farm scenario" begin
#             # Simulate a more realistic wind farm control scenario
#             # 5 time steps, 4 turbines with varying induction factors
#             con_induction_data = [
#                 0.0   0.30  0.28  0.32  0.29;  # t=0s
#                 10.0  0.31  0.29  0.33  0.30;  # t=10s  
#                 20.0  0.33  0.31  0.35  0.32;  # t=20s
#                 30.0  0.32  0.30  0.34  0.31;  # t=30s
#                 40.0  0.30  0.28  0.32  0.29   # t=40s
#             ]
            
#             # Test interpolation at intermediate times
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 1, 5.0) ≈ 0.305
#             @test getInduction(FLORIDyn.Induction_MPC(), con_induction_data, 2, 15.0) ≈ 0.30
            
#             # Test multiple turbines at once
#             result = getInduction(FLORIDyn.Induction_MPC(), con_induction_data, [1, 2, 3, 4], 25.0)
#             expected = [0.325, 0.305, 0.345, 0.315]
#             @test result ≈ expected
#         end
        
#         @testset "Performance and numerical stability" begin
#             # Test with larger dataset to ensure performance
#             n_times = 100
#             n_turbines = 10
            
#             # Generate test data
#             times = collect(0.0:1.0:(n_times-1))
#             induction_data = zeros(n_times, n_turbines + 1)
#             induction_data[:, 1] = times
            
#             for i in 1:n_turbines
#                 # Create varying induction factors for each turbine
#                 base_induction = 0.25 + 0.1 * (i-1) / (n_turbines-1)
#                 induction_data[:, i+1] = base_induction .+ 0.05 * sin.(2π * times / 50.0)
#             end
            
#             # Test performance with large dataset
#             @test typeof(getInduction(FLORIDyn.Induction_MPC(), induction_data, 1, 25.5)) == Float64
#             @test typeof(getInduction(FLORIDyn.Induction_MPC(), induction_data, 1:5, 25.5)) == Vector{Float64}
            
#             # Test numerical stability near boundaries
#             @test isfinite(getInduction(FLORIDyn.Induction_MPC(), induction_data, 1, 0.0))
#             @test isfinite(getInduction(FLORIDyn.Induction_MPC(), induction_data, 1, n_times-1))
#         end
#     end
    
#     @testset "Type consistency between getYaw and getInduction" begin
#         # Ensure both functions have consistent behavior and return types
#         yaw_data = [0.0 10.0; 1.0 15.0; 2.0 20.0]
#         induction_data = [0.0 0.30; 1.0 0.32; 2.0 0.34]
        
#         # Test single turbine returns
#         yaw_result = getYaw(FLORIDyn.Yaw_SOWFA(), yaw_data, 1, 1.5)
#         induction_result = getInduction(FLORIDyn.Induction_MPC(), induction_data, 1, 1.5)
#         @test typeof(yaw_result) == typeof(induction_result) == Float64
        
#         # Test multiple turbines returns
#         yaw_multi = getYaw(FLORIDyn.Yaw_SOWFA(), yaw_data, [1], 1.5)
#         induction_multi = getInduction(FLORIDyn.Induction_MPC(), induction_data, [1], 1.5)
#         @test typeof(yaw_multi) == typeof(induction_multi) == Vector{Float64}
#     end
end

nothing
