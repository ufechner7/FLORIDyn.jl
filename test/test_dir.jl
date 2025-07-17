# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using Random

Random.seed!(1234)

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
        phi = getWindDirT(dir_mode, WindDir, iT)
        result = [269.0527533560426, 270.54014974707036, 269.78339785902375]
        for (i, ph) in pairs(phi)
            @test ph ≈ result[i]
        end
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

        # Call the function
        phi = getWindDirT(dir_mode, WindDir, iT, t)
        @test size(phi) == (2,1)
        @test phi[1] ≈ 25.84752589085377
        @test phi[2] ≈ 25.140544918823198128
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
        phi = getWindDirT(dir_mode, WindDir, 1, 12.5)

        @test length(phi) == 1
        @test phi[1] ≈ 355.84437113031197
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

        # Call the function
        phi = getWindDirT(dir_mode, WindDirNow, WindDir)
        @test size(phi) == (3,1)
        @test phi[1] ≈ 9.80509368889485
        @test phi[2] ≈ 21.48940455337165
        @test phi[3] ≈ 31.437827702779927
    end
end

nothing
