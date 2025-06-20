# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using Random

Random.seed!(1234)

struct WindDirType
    Data::Float64
    CholSig::Matrix{Float64}
end

# Define a struct for WindDir
struct WindDirMatrix
    Data::Matrix{Float64}      # Nx2 matrix: column 1 = time, column 2 = phi
    CholSig::Matrix{Float64}   # Cholesky factor of covariance matrix (nT x nT)
end

struct WindDirTriple
    Init::Vector{Float64}
    CholSig::Matrix{Float64}
    MeanPull::Float64
end

@testset "windfield.jl" begin
    dir_mode = Direction_Constant()
    WindDir = 270
    iT = [1, 2, 3]
    phi = getWindDirT(dir_mode, WindDir, iT, nothing)
    for ph in phi
        @test ph ≈ 270.0
    end

    #########################################

    dir_mode = Direction_Constant_wErrorCov()
    WindDir = WindDirType(270.0, cholesky(Matrix{Float64}(I, 3, 3)).L)
    iT = [1, 2, 3]
    phi = getWindDirT(dir_mode, WindDir, iT)
    result = [269.0527533560426, 270.54014974707036, 269.78339785902375]
    for (i, ph) in pairs(phi)
        @test ph ≈ result[i]
    end

    #########################################

    dir_mode = Direction_EnKF_InterpTurbine()

    # Suppose WindDir is a matrix where each row is [time, phi_T0, phi_T1, ...]
    WindDir = [
        0.0  10.0  20.0
        1.0  12.0  22.0
        2.0  14.0  24.0
    ]
    phi = getWindDirT_EnKF(dir_mode, WindDir, 1, 0.5)
    @test phi ≈ 11.0

    dir_mode = Direction_Interpolation()
    phi = getWindDirT(dir_mode, WindDir, 1, 0.5)
    @test phi[1] ≈ 11.0

    ##############################################

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

    ###################################

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

    ##############################################

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
    phi = getWindDirT(Direction_InterpTurbine_wErrorCov(), WindDir, 1, 12.5)

    @test length(phi) == 1
    @test phi[1] ≈ 355.84437113031197

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
    phi = getWindDirT(Direction_RW_with_Mean(), WindDirNow, WindDir)

end
nothing
