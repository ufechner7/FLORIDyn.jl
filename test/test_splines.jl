# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test

# Import the function to test
import FLORIDyn: interpolate_hermite_spline

@testset verbose=true "Hermite Spline Interpolation Tests" begin
    
    @testset "Equidistant control points - Basic functionality" begin
        @testset "Two-point linear interpolation" begin
            correction = [0.0, 1.0]
            
            # Test endpoints
            @test interpolate_hermite_spline(0.0, correction) ≈ 0.0
            @test interpolate_hermite_spline(1.0, correction) ≈ 1.0
            
            # Test midpoint (should be linear)
            @test interpolate_hermite_spline(0.5, correction) ≈ 0.5
            
            # Test quarter points
            @test interpolate_hermite_spline(0.25, correction) ≈ 0.25
            @test interpolate_hermite_spline(0.75, correction) ≈ 0.75
        end
        
        @testset "Three-point interpolation" begin
            correction = [0.0, 1.0, 0.0]
            
            # Test endpoints
            @test interpolate_hermite_spline(0.0, correction) ≈ 0.0
            @test interpolate_hermite_spline(1.0, correction) ≈ 0.0
            
            # Test middle control point
            @test interpolate_hermite_spline(0.5, correction) ≈ 1.0
            
            # Due to Hermite spline, values between control points should be smooth
            # At quarter point, should be between 0.0 and 1.0
            val_quarter = interpolate_hermite_spline(0.25, correction)
            @test 0.0 < val_quarter < 1.0
            
            # At three-quarter point, should be between 0.0 and 1.0
            val_three_quarter = interpolate_hermite_spline(0.75, correction)
            @test 0.0 < val_three_quarter < 1.0
        end
        
        @testset "Four-point interpolation" begin
            correction = [0.0, 1.0, 2.0, 1.0]
            
            # Test all control points
            @test interpolate_hermite_spline(0.0, correction) ≈ 0.0
            @test interpolate_hermite_spline(1.0/3.0, correction) ≈ 1.0
            @test interpolate_hermite_spline(2.0/3.0, correction) ≈ 2.0
            @test interpolate_hermite_spline(1.0, correction) ≈ 1.0
            
            # Test intermediate values
            val_mid_first = interpolate_hermite_spline(1.0/6.0, correction)
            @test 0.0 < val_mid_first < 1.0
        end
        
        @testset "Constant values" begin
            correction = [5.0, 5.0, 5.0, 5.0]
            
            # All points should return the constant value
            @test interpolate_hermite_spline(0.0, correction) ≈ 5.0
            @test interpolate_hermite_spline(0.25, correction) ≈ 5.0 atol=1e-10
            @test interpolate_hermite_spline(0.5, correction) ≈ 5.0 atol=1e-10
            @test interpolate_hermite_spline(0.75, correction) ≈ 5.0 atol=1e-10
            @test interpolate_hermite_spline(1.0, correction) ≈ 5.0
        end
        
        @testset "Monotonicity preservation" begin
            # Increasing values
            correction = [0.0, 1.0, 2.0, 3.0]
            vals = [interpolate_hermite_spline(s, correction) for s in 0.0:0.1:1.0]
            
            # Check that values are generally increasing (allowing for small numerical variations)
            for i in 2:length(vals)
                @test vals[i] >= vals[i-1] - 1e-10
            end
        end
    end
    
    @testset "Non-equidistant control points - Basic functionality" begin
        @testset "Two-point linear interpolation" begin
            correction = [0.0, 1.0]
            s_positions = [0.0, 1.0]
            
            # Test endpoints
            @test interpolate_hermite_spline(0.0, correction, s_positions) ≈ 0.0
            @test interpolate_hermite_spline(1.0, correction, s_positions) ≈ 1.0
            
            # Test midpoint
            @test interpolate_hermite_spline(0.5, correction, s_positions) ≈ 0.5
        end
        
        @testset "Three-point non-equidistant" begin
            correction = [0.0, 1.0, 0.0]
            s_positions = [0.0, 0.3, 1.0]  # Non-equidistant
            
            # Test control points
            @test interpolate_hermite_spline(0.0, correction, s_positions) ≈ 0.0
            @test interpolate_hermite_spline(0.3, correction, s_positions) ≈ 1.0
            @test interpolate_hermite_spline(1.0, correction, s_positions) ≈ 0.0
            
            # Test intermediate values
            val_before_peak = interpolate_hermite_spline(0.15, correction, s_positions)
            @test 0.0 < val_before_peak < 1.0
            
            val_after_peak = interpolate_hermite_spline(0.65, correction, s_positions)
            @test 0.0 < val_after_peak < 1.0
        end
        
        @testset "Four-point non-equidistant" begin
            correction = [0.0, 2.0, 1.0, 3.0]
            s_positions = [0.0, 0.2, 0.5, 1.0]
            
            # Test control points
            @test interpolate_hermite_spline(0.0, correction, s_positions) ≈ 0.0
            @test interpolate_hermite_spline(0.2, correction, s_positions) ≈ 2.0
            @test interpolate_hermite_spline(0.5, correction, s_positions) ≈ 1.0
            @test interpolate_hermite_spline(1.0, correction, s_positions) ≈ 3.0
        end
        
        @testset "Heavily skewed spacing" begin
            correction = [0.0, 1.0, 2.0]
            s_positions = [0.0, 0.9, 1.0]  # Heavily skewed to the right
            
            # Test control points
            @test interpolate_hermite_spline(0.0, correction, s_positions) ≈ 0.0
            @test interpolate_hermite_spline(0.9, correction, s_positions) ≈ 1.0
            @test interpolate_hermite_spline(1.0, correction, s_positions) ≈ 2.0
            
            # Most of the curve should have values close to 0
            val_mid = interpolate_hermite_spline(0.45, correction, s_positions)
            @test val_mid < 1.0
        end
    end
    
    @testset "Boundary conditions and edge cases" begin
        @testset "Values at s=0 and s=1" begin
            correction = [1.5, 2.5, 3.5, 4.5]
            
            @test interpolate_hermite_spline(0.0, correction) ≈ 1.5
            @test interpolate_hermite_spline(1.0, correction) ≈ 4.5
        end
        
        @testset "Values slightly outside [0,1] should be clamped" begin
            correction = [0.0, 1.0, 2.0]
            
            # Values at exactly 0 and 1
            @test interpolate_hermite_spline(0.0, correction) ≈ 0.0
            @test interpolate_hermite_spline(1.0, correction) ≈ 2.0
            
            # Test that values very close to boundaries work
            @test interpolate_hermite_spline(0.0001, correction) >= 0.0
            @test interpolate_hermite_spline(0.9999, correction) <= 2.0 + 0.1
        end
        
        @testset "Negative values" begin
            correction = [-1.0, -2.0, -3.0]
            
            @test interpolate_hermite_spline(0.0, correction) ≈ -1.0
            @test interpolate_hermite_spline(0.5, correction) ≈ -2.0
            @test interpolate_hermite_spline(1.0, correction) ≈ -3.0
        end
    end
    
    @testset "Smoothness properties (C1 continuity)" begin
        @testset "Continuous first derivative at control points" begin
            correction = [0.0, 1.0, 0.5, 2.0]
            
            # Numerically check derivative continuity at control point
            ε = 1e-6
            s_test = 1.0/3.0  # Second control point
            
            # Left derivative
            val_left = interpolate_hermite_spline(s_test - ε, correction)
            val_center = interpolate_hermite_spline(s_test, correction)
            deriv_left = (val_center - val_left) / ε
            
            # Right derivative
            val_right = interpolate_hermite_spline(s_test + ε, correction)
            deriv_right = (val_right - val_center) / ε
            
            # Derivatives should be approximately equal (C1 continuity)
            @test deriv_left ≈ deriv_right rtol=0.01
        end
    end
    
    @testset "Assertion tests for invalid inputs" begin
        @testset "Equidistant version - insufficient points" begin
            @test_throws AssertionError interpolate_hermite_spline(0.5, [1.0])
            @test_throws AssertionError interpolate_hermite_spline(0.5, Float64[])
        end
        
        @testset "Non-equidistant version - mismatched lengths" begin
            correction = [0.0, 1.0, 2.0]
            s_positions = [0.0, 0.5]  # Wrong length
            
            @test_throws AssertionError interpolate_hermite_spline(0.5, correction, s_positions)
        end
        
        @testset "Non-equidistant version - invalid position range" begin
            correction = [0.0, 1.0, 2.0]
            s_positions = [0.1, 0.5, 1.0]  # Doesn't start at 0
            
            @test_throws AssertionError interpolate_hermite_spline(0.5, correction, s_positions)
            
            s_positions = [0.0, 0.5, 0.9]  # Doesn't end at 1
            @test_throws AssertionError interpolate_hermite_spline(0.5, correction, s_positions)
        end
        
        @testset "Non-equidistant version - insufficient points" begin
            correction = [1.0]
            s_positions = [0.5]
            
            @test_throws AssertionError interpolate_hermite_spline(0.5, correction, s_positions)
        end
    end
    
    @testset "Comparison between equidistant and non-equidistant versions" begin
        @testset "Should match when positions are equidistant" begin
            correction = [0.0, 1.0, 2.0, 1.5]
            n_points = length(correction)
            s_positions = [i/(n_points-1) for i in 0:(n_points-1)]
            
            # Test at multiple points
            for s in 0.0:0.1:1.0
                val_eq = interpolate_hermite_spline(s, correction)
                val_noneq = interpolate_hermite_spline(s, correction, s_positions)
                @test val_eq ≈ val_noneq atol=1e-10
            end
        end
    end
    
    @testset "Stress tests with many control points" begin
        @testset "Ten control points" begin
            correction = sin.(range(0, 2π, length=10))
            
            # Test all control points
            for i in 0:9
                s = i / 9.0
                val = interpolate_hermite_spline(s, correction)
                @test val ≈ correction[i+1] atol=1e-10
            end
            
            # Test intermediate values are bounded
            for s in 0.0:0.01:1.0
                val = interpolate_hermite_spline(s, correction)
                @test -1.5 <= val <= 1.5
            end
        end
        
        @testset "Twenty control points non-equidistant" begin
            n_points = 20
            correction = rand(n_points)
            s_positions = sort([0.0; cumsum(rand(n_points-2)); 1.0])
            s_positions ./= s_positions[end]  # Normalize to [0,1]
            
            # Test all control points
            for i in 1:n_points
                val = interpolate_hermite_spline(s_positions[i], correction, s_positions)
                @test val ≈ correction[i] atol=1e-10
            end
        end
    end
end
