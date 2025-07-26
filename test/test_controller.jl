# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, Interpolations

# Import the getYaw function directly
import FLORIDyn: getYaw

@testset "getYaw function tests" begin
    
    @testset "Basic interpolation with single turbine" begin
        # Test data: time column + one turbine yaw column
        ConYawData = [
            0.0  10.0;
            1.0  15.0;
            2.0  20.0;
            3.0  25.0;
            4.0  30.0
        ]
        
        # Test interpolation at exact time points
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 0.0) ≈ 10.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 1.0) ≈ 15.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 2.0) ≈ 20.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 3.0) ≈ 25.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 4.0) ≈ 30.0
        
        # Test interpolation at intermediate time points
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 0.5) ≈ 12.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 1.5) ≈ 17.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 2.5) ≈ 22.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 3.5) ≈ 27.5
    end
    
    @testset "Multiple turbines interpolation" begin
        # Test data: time column + three turbine yaw columns
        ConYawData = [
            0.0  10.0  5.0   0.0;
            1.0  15.0  10.0  5.0;
            2.0  20.0  15.0  10.0;
            3.0  25.0  20.0  15.0
        ]
        
        # Test single turbine access
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 1.5) ≈ 17.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 2, 1.5) ≈ 12.5
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 3, 1.5) ≈ 7.5
        
        # Test multiple turbine access
        result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [1, 2, 3], 1.5)
        @test result ≈ [17.5, 12.5, 7.5]
        
        # Test subset of turbines
        result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [1, 3], 2.5)
        @test result ≈ [22.5, 12.5]
        
        # Test single turbine in vector format
        result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [2], 0.5)
        @test result ≈ [7.5]
    end
    
    @testset "Out of bounds time handling" begin
        ConYawData = [
            1.0  10.0;
            2.0  20.0;
            3.0  30.0;
            4.0  40.0
        ]
        
        # Test time before first time point (should clamp to first time and issue warning)
        @test_logs (:warn, "The time 0.5 is out of bounds, will use 1.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 0.5)
            @test result ≈ 10.0
        end
        
        # Test time after last time point (should clamp to last time and issue warning)
        @test_logs (:warn, "The time 5.0 is out of bounds, will use 4.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 5.0)
            @test result ≈ 40.0
        end
        
        # Test with multiple turbines and out of bounds time
        ConYawData_multi = [
            1.0  10.0  5.0;
            2.0  20.0  15.0;
            3.0  30.0  25.0
        ]
        
        @test_logs (:warn, "The time 0.0 is out of bounds, will use 1.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData_multi, [1, 2], 0.0)
            @test result ≈ [10.0, 5.0]
        end
    end
    
    @testset "Edge cases and error handling" begin
        ConYawData = [
            0.0  10.0  5.0;
            1.0  20.0  15.0;
            2.0  30.0  25.0
        ]
        
        # Test invalid turbine index types
        @test_throws ErrorException getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1.5, 1.0)
        @test_throws ErrorException getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, "1", 1.0)
        @test_throws ErrorException getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [1.5, 2.0], 1.0)
        
        # Test valid edge cases
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, Int32(1), 1.0) ≈ 20.0  # Different integer type
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [Int32(1), Int64(2)], 1.0) ≈ [20.0, 15.0]  # Mixed integer types
    end
    
    @testset "Single time point data" begin
        # Test with only one time point
        ConYawData = [5.0  25.0  35.0]
        
        # Should return the single value regardless of requested time
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 5.0) ≈ 25.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 2, 5.0) ≈ 35.0
        
        # Out of bounds should still issue warnings but return the single value
        @test_logs (:warn, "The time 0.0 is out of bounds, will use 5.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 0.0)
            @test result ≈ 25.0
        end
        
        @test_logs (:warn, "The time 10.0 is out of bounds, will use 5.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [1, 2], 10.0)
            @test result ≈ [25.0, 35.0]
        end
    end
    
    @testset "Two time point data" begin
        # Test with exactly two time points (minimum for linear interpolation)
        ConYawData = [
            0.0  0.0   10.0;
            10.0 100.0 20.0
        ]
        
        # Test exact endpoints
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 0.0) ≈ 0.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 10.0) ≈ 100.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 2, 0.0) ≈ 10.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 2, 10.0) ≈ 20.0
        
        # Test interpolation
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 5.0) ≈ 50.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 2, 5.0) ≈ 15.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [1, 2], 2.5) ≈ [25.0, 12.5]
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
        
        # Test at specific time points
        result_t150 = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1:9, 150.0)
        expected_t150 = [2.5, 7.5, 12.5, 7.5, 12.5, 17.5, 12.5, 17.5, 22.5]
        @test result_t150 ≈ expected_t150
        
        # Test subset of turbines at different time
        result_subset = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [1, 5, 9], 250.0)
        expected_subset = [7.5, 17.5, 27.5]
        @test result_subset ≈ expected_subset
        
        # Test single turbine across different times
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 5, 0.0) ≈ 5.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 5, 300.0) ≈ 20.0
        @test getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 5, 600.0) ≈ 0.0
    end
    
    @testset "Extrapolation behavior (Flat boundary condition)" begin
        # The function uses Flat() extrapolation, so values beyond bounds should be constant
        ConYawData = [
            1.0  10.0;
            2.0  20.0;
            3.0  30.0
        ]
        
        # Note: The function clamps time to bounds and issues warnings,
        # but the underlying interpolation would use flat extrapolation
        # We test this by checking the clamping behavior
        
        @test_logs (:warn, "The time 0.5 is out of bounds, will use 1.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 0.5)
            @test result ≈ 10.0  # Should be the first value
        end
        
        @test_logs (:warn, "The time 4.0 is out of bounds, will use 3.0 instead.") begin
            result = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1, 4.0)
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
        
        # Test single turbine access
        @time result1 = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 25, 500.5)
        @test isa(result1, Float64)
        
        # Test multiple turbine access
        @time result_all = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, 1:n_turbines, 500.5)
        @test length(result_all) == n_turbines
        @test all(isa.(result_all, Float64))
        
        # Test subset access
        @time result_subset = getYaw(FLORIDyn.Yaw_SOWFA(), ConYawData, [1, 10, 25, 40, 50], 750.3)
        @test length(result_subset) == 5
    end
end
