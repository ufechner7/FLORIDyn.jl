# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using Random
using Statistics

FLORIDyn.set_rng(MersenneTwister(1234))

@testset verbose = true "wind dir                                                " begin
    @testset "getWindDirT(Direction_Constant(), ...)" begin
        dir_mode = Direction_Constant()
        WindDir = 270.0
        iT = [1, 2, 3]
        phi = getWindDirT(dir_mode, WindDir, iT, nothing)
        for ph in phi
            @test ph ≈ 270.0
        end
    end

    @testset "getWindDirT(Direction_Constant_wErrorCov(), ...)" begin
        dir_mode = Direction_Constant_wErrorCov()
        WindDir = WindDirType(270.0, cholesky(Matrix{Float64}(I, 3, 3)).L)
        iT = [1, 2, 3]
        
        # Test single call
        phi = getWindDirT(dir_mode, WindDir, iT, 0.0)
        @test length(phi) == 3
        @test all(isfinite.(phi))
        
        # Test statistical properties with multiple samples
        n_samples = 1000
        samples = [getWindDirT(dir_mode, WindDir, iT, 0.0) for _ in 1:n_samples]
        sample_matrix = hcat(samples...)  # 3 x n_samples matrix
        
        # Test that means are close to the base wind direction
        means = mean(sample_matrix, dims=2)
        @test all(abs.(means .- 270.0) .< 0.5)  # Should be close to 270.0
        
        # Test that variances are reasonable (identity covariance should give ~1.0 variance)
        vars = var(sample_matrix, dims=2)
        @test all(0.5 .< vars .< 2.0)  # Should be around 1.0
    end

    @testset "getWindDirT_EnKF(Direction_EnKF_InterpTurbine(), ...)" begin
        dir_mode = Direction_EnKF_InterpTurbine()

        # Suppose WindDir is a matrix where each row is [time, phi_T0, phi_T1, ...]
        WindDir = [
            0.0  10.0  20.0
            1.0  12.0  22.0
            2.0  14.0  24.0
        ]
        phi = getWindDirT_EnKF(dir_mode, WindDir, 1, 0.5)
        @test phi ≈ 11.0
    end

    @testset "getWindDirT_EnKF(Direction_Interpolation(), ...)" begin
        dir_mode = Direction_Interpolation()
        phi = getWindDirT(dir_mode, WindDir, 1, 0.5)
        @test phi[1] ≈ 11.0
    end

    # Example wind direction data: times in first column, turbines in next columns
    WindDir = [
        0.0  10.0  20.0;
        1.0  15.0  25.0;
        2.0  20.0  30.0
    ]

    @testset "getWindDirT(Direction_InterpTurbine(), ...)" begin
        dir_mode = Direction_InterpTurbine()
        # Test for interpolation at t = 0.5 (between 0.0 and 1.0)
        @test getWindDirT(dir_mode, WindDir, 1, 0.5) ≈ 12.5
        @test getWindDirT(dir_mode, WindDir, 2, 0.5) ≈ 22.5

        # Test for exact time at t = 1.0
        @test getWindDirT(dir_mode, WindDir, 1, 1.0) == 15.0
        @test getWindDirT(dir_mode, WindDir, 2, 1.0) == 25.0

        # Test for time below minimum (should clamp to tmin)
        @test getWindDirT(dir_mode, WindDir, 1, -1.0) == 10.0

        # Test for time above maximum (should clamp to tmax)
        @test getWindDirT(dir_mode, WindDir, 2, 5.0) == 30.0
    end

    @testset "getWindDirT(Direction_Interpolation_wErrorCov(), ...)" begin
        dir_mode = Direction_Interpolation_wErrorCov()

        # Example wind direction data (time, phi)
        wind_data = [
            0.0  10.0;
            5.0  20.0;
            10.0 30.0
        ]

        # Example Cholesky factor (for 2 turbines)
        chol_sig = cholesky([1.0 0.5; 0.5 1.0]).L

        # Create WindDir instance
        WindDir = WindDirMatrix(wind_data, chol_sig)

        # Example turbine indices (for 2 turbines)
        iT = [1, 2]

        # Example time
        t = 7.0

        # Call the function and test statistical properties
        n_samples = 500
        samples = [getWindDirT(dir_mode, WindDir, iT, t) for _ in 1:n_samples]
        sample_matrix = hcat(samples...)  # 2 x n_samples matrix
        
        # Test basic properties
        phi = getWindDirT(dir_mode, WindDir, iT, t)
        @test size(phi) == (2,1)  # Corrected expected size
        @test all(isfinite.(phi))
        
        # Test that the mean is close to the interpolated base value (26.0 at t=7.0)
        expected_base = 10.0 + (20.0 - 10.0) * (7.0 - 0.0) / (5.0 - 0.0)  # ≈ 26.0
        means = mean(sample_matrix, dims=2)
        @test all(abs.(means .- expected_base) .< 1.0)
        
        # Test that variances are reasonable
        vars = var(sample_matrix, dims=2)
        @test all(0.5 .< vars .< 3.0)
    end

    @testset "getWindDirT(Direction_InterpTurbine(), ...)" begin
        dir_mode = Direction_InterpTurbine()
        # Example wind direction data:
        # Columns: time, phi_T0, phi_T1, phi_T2
        WindDir = [
            0.0   10.0  20.0  30.0;
            10.0  15.0  25.0  35.0;
            20.0  20.0  30.0  40.0
        ]

        # Turbine index (1-based): e.g., 2 for phi_T1
        iT = 2

        # Time at which to interpolate wind direction
        t = 5.0

        # Call the function
        phi = getWindDirT(dir_mode, WindDir, iT, t)
        @test phi ≈ 22.5
    end

    @testset "getWindDirT(Direction_InterpTurbine_wErrorCov(), ...)" begin

        dir_mode = Direction_InterpTurbine_wErrorCov()
        
        # Example wind direction data (time, phi)
        wind_data = [
            0.0  10.0;
            5.0  20.0;
            10.0 30.0
        ]

        # Example Cholesky factor (for 2 turbines)
        chol_sig = cholesky([1.0 0.5; 0.5 1.0]).L

        # Create WindDir instance
        WindDir = WindDirMatrix(wind_data, chol_sig)

        # Example time series: 0, 10, 20 seconds
        times = [0.0, 10.0, 20.0]
        # Example wind directions for 2 turbines at each time
        phi_T0 = [350.0, 355.0, 360.0]
        phi_T1 = [10.0, 15.0, 20.0]
        # Combine into Data matrix: each row is [time, phi_T0, phi_T1]
        Data = hcat(times, phi_T0, phi_T1)

        # Example covariance matrix and its Cholesky factor
        cov = [1.0 0.5; 0.5 1.0]
        CholSig = cholesky(cov).L

        # Create WindDir struct
        WindDir = WindDirMatrix(Data, CholSig)
        
        # Test single call
        phi = getWindDirT(dir_mode, WindDir, 1, 12.5)
        @test length(phi) == 1
        @test isfinite(phi[1])
        
        # Test statistical properties with multiple samples
        n_samples = 500
        samples = [getWindDirT(dir_mode, WindDir, 1, 12.5)[1] for _ in 1:n_samples]
        
        # Expected interpolated value at t=12.5 for turbine 1: 355 + (360-355)*(12.5-10)/(20-10) = 356.25
        expected_base = 355.0 + (360.0 - 355.0) * (12.5 - 10.0) / (20.0 - 10.0)
        sample_mean = mean(samples)
        sample_var = var(samples)
        
        # Test that mean is close to expected interpolated value
        @test abs(sample_mean - expected_base) < 0.5
        
        # Test that variance is reasonable
        @test 0.5 < sample_var < 3.0 
    end

    @testset "getWindDirT(Direction_RW_with_Mean(), ...)" begin
        dir_mode = Direction_RW_with_Mean()

        # Suppose we have 3 turbines
        WindDirNow = [10.0, 20.0, 30.0]           # Current wind directions (degrees)
        Init = [15.0, 25.0, 35.0]                 # Mean wind directions (degrees)
        CholSig = [1.0 0.2 0.1;                   # Cholesky factor of covariance matrix
                0.0 1.0 0.3;
                0.0 0.0 1.0]
        MeanPull = 0.05                           # Mean reversion factor

        # Create WindDir struct
        WindDir = WindDirTriple(Init, CholSig, MeanPull)

        # Test single call
        phi = getWindDirT(dir_mode, WindDirNow, WindDir)
        @test size(phi) == (3,1)
        @test all(isfinite.(phi))
        
        # Test statistical properties with multiple samples
        n_samples = 500
        samples = [getWindDirT(dir_mode, WindDirNow, WindDir) for _ in 1:n_samples]
        sample_matrix = hcat(samples...)  # 3 x n_samples matrix
        
        # Expected direction after one step: WindDirNow + MeanPull * (Init - WindDirNow)
        # For the random walk with mean reversion
        expected_directions = WindDirNow .+ MeanPull .* (Init .- WindDirNow)
        
        means = mean(sample_matrix, dims=2)
        
        # Test that means are close to expected directions (accounting for random noise)
        @test all(abs.(means .- expected_directions) .< 1.0)
        
        # Test that variances are reasonable
        vars = var(sample_matrix, dims=2)
        @test all(0.2 .< vars .< 5.0)
    end
end

nothing
