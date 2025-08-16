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
        buffer = Vector{Float64}(undef, 1)
        phi = getWindDirT!(buffer, dir_mode, WindDir, 1, 0.5)
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

    @testset "getWindDirT(Direction_RW_with_Mean(), WindDirTriple, iT, t)" begin
        dir_mode = Direction_RW_with_Mean()

        # Create test data
        Init = [15.0, 25.0, 35.0, 45.0]          # Mean wind directions for 4 turbines (degrees)
        CholSig = [1.0 0.2 0.1 0.05;             # Cholesky factor of covariance matrix (4x4)
                   0.0 1.0 0.3 0.15;
                   0.0 0.0 1.0 0.25;
                   0.0 0.0 0.0 1.0]
        MeanPull = 0.1                           # Mean reversion factor

        # Create WindDir struct
        WindDir = WindDirTriple(Init, CholSig, MeanPull)

        @testset "Single turbine selection" begin
            # Test with single turbine index
            iT = 2
            t = 5.0  # Time value (unused in current implementation)
            
            # Test single call
            phi = getWindDirT(dir_mode, WindDir, iT, t)
            @test isa(phi, Float64)  # Should return scalar for single turbine
            @test isfinite(phi)
            
            # Test statistical properties with multiple samples
            n_samples = 1000
            samples = [getWindDirT(dir_mode, WindDir, iT, t) for _ in 1:n_samples]
            
            # Expected direction: Init[iT] + MeanPull * (Init[iT] - Init[iT]) = Init[iT]
            # Since WindDirNow = Init[iT] initially
            expected_mean = Init[iT]
            sample_mean = mean(samples)
            sample_var = var(samples)
            
            # Test that mean is close to initial value (no mean reversion when starting at mean)
            @test abs(sample_mean - expected_mean) < 0.3
            
            # Test that variance is reasonable (should be around CholSig[iT,iT]^2)
            expected_var = CholSig[iT, iT]^2
            @test 0.5 * expected_var < sample_var < 2.0 * expected_var
        end

        @testset "Multiple turbine selection" begin
            # Test with multiple turbine indices
            iT = [1, 3, 4]
            t = 10.0  # Time value (unused in current implementation)
            
            # Test single call
            phi = getWindDirT(dir_mode, WindDir, iT, t)
            @test isa(phi, Matrix{Float64})
            @test size(phi) == (3, 1)  # Column vector format
            @test all(isfinite.(phi))
            
            # Test statistical properties with multiple samples
            n_samples = 800
            samples = [getWindDirT(dir_mode, WindDir, iT, t) for _ in 1:n_samples]
            sample_matrix = hcat([vec(s) for s in samples]...)  # Convert each column vector to regular vector, then combine
            
            # Expected directions: Init[iT] for each selected turbine
            expected_means = Init[iT]
            sample_means = mean(sample_matrix, dims=2)[:, 1]
            
            # Test that means are close to initial values
            @test all(abs.(sample_means .- expected_means) .< 0.3)
            
            # Test that variances are reasonable
            sample_vars = var(sample_matrix, dims=2)[:, 1]
            expected_vars = [CholSig[i, i]^2 for i in iT]
            @test all(0.5 .* expected_vars .< sample_vars .< 2.0 .* expected_vars)
        end

        @testset "All turbines selection" begin
            # Test with all turbine indices
            iT = [1, 2, 3, 4]
            t = 15.0
            
            # Test single call
            phi = getWindDirT(dir_mode, WindDir, iT, t)
            @test isa(phi, Matrix{Float64})
            @test size(phi) == (4, 1)  # Column vector format
            @test all(isfinite.(phi))
            
            # Test statistical properties
            n_samples = 500
            samples = [getWindDirT(dir_mode, WindDir, iT, t) for _ in 1:n_samples]
            sample_matrix = hcat([vec(s) for s in samples]...)  # Convert each column vector to regular vector, then combine
            
            # Test means
            expected_means = Init
            sample_means = mean(sample_matrix, dims=2)[:, 1]
            @test all(abs.(sample_means .- expected_means) .< 0.3)
            
            # Test that covariance structure is preserved
            sample_cov = cov(sample_matrix, dims=2)
            # The sample covariance should be approximately CholSig * CholSig'
            expected_cov = CholSig * CholSig'
            
            # Test diagonal elements (variances)
            for i in 1:4
                @test 0.5 * expected_cov[i, i] < sample_cov[i, i] < 2.0 * expected_cov[i, i]
            end
            
            # Test that off-diagonal correlations have the right sign
            @test sign(sample_cov[1, 2]) == sign(expected_cov[1, 2])
            @test sign(sample_cov[2, 3]) == sign(expected_cov[2, 3])
        end

        @testset "Edge cases" begin
            # Test with zero mean reversion
            WindDir_zero = WindDirTriple(Init, CholSig, 0.0)
            iT = 2
            t = 0.0
            
            phi = getWindDirT(dir_mode, WindDir_zero, iT, t)
            @test isfinite(phi)
            
            # Test with high mean reversion
            WindDir_high = WindDirTriple(Init, CholSig, 1.0)
            phi_high = getWindDirT(dir_mode, WindDir_high, iT, t)
            @test isfinite(phi_high)
            
            # Test with single turbine (iT = 1)
            iT_single = 1
            phi_single = getWindDirT(dir_mode, WindDir, iT_single, t)
            @test isa(phi_single, Float64)
            @test isfinite(phi_single)
        end

        @testset "Time independence" begin
            # The current implementation doesn't use time `t`, so results should be identical for different times
            # (though random, so we test the distribution properties)
            iT = [2, 4]
            t1 = 0.0
            t2 = 100.0
            
            n_samples = 300
            samples_t1 = [getWindDirT(dir_mode, WindDir, iT, t1) for _ in 1:n_samples]
            samples_t2 = [getWindDirT(dir_mode, WindDir, iT, t2) for _ in 1:n_samples]
            
            # Convert to regular matrices for statistical analysis
            matrix_t1 = hcat([vec(s) for s in samples_t1]...)
            matrix_t2 = hcat([vec(s) for s in samples_t2]...)
            
            # Means should be similar (both should be around Init[iT])
            mean_t1 = mean(matrix_t1, dims=2)[:, 1]
            mean_t2 = mean(matrix_t2, dims=2)[:, 1]
            expected = Init[iT]
            
            @test all(abs.(mean_t1 .- expected) .< 0.4)
            @test all(abs.(mean_t2 .- expected) .< 0.4)
            
            # Variances should be similar
            var_t1 = var(matrix_t1, dims=2)[:, 1]
            var_t2 = var(matrix_t2, dims=2)[:, 1]
            
            # They should be of similar magnitude (within factor of 2)
            @test all(var_t1 ./ var_t2 .< 2.0)
            @test all(var_t2 ./ var_t1 .< 2.0)
        end
    end
end

nothing
