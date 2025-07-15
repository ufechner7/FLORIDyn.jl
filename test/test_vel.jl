# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using Statistics
using Random

@testset "wind velocity  " begin
    @testset "getWindSpeedT(Velocity_Constant(), ...)" begin
        vc = Velocity_Constant()

        @testset "Scalar index test" begin
            result = getWindSpeedT(vc, 10.0, 1)
            @test result ≈ 10.0
            @test typeof(result) == Float64
        end

        @testset "Array of indices test" begin
            result = getWindSpeedT(vc, 5.0, [1,2,3])
            @test result ≈ [5.0, 5.0, 5.0]
            @test typeof(result) == Vector{Float64}
        end

        @testset "Higher-dimensional index array" begin
            idx = [1 2; 3 4]
            result = getWindSpeedT(vc, 7.5, idx)
            @test result ≈ [7.5 7.5; 7.5 7.5]
            @test size(result) == (2, 2)
        end

        @testset "WindVel as Int" begin
            result = getWindSpeedT(vc, 3, [10, 20])
            @test result == [3, 3]
            @test eltype(result) == Int
        end

        @testset "Empty index array" begin
            result = getWindSpeedT(vc, 6.0, Int[])
            @test result == Float64[]
            @test length(result) == 0
        end
    end
    @testset "get_wind_speed_t Unit Tests" begin
        # Set fixed random seed for consistent results
        Random.seed!(1234)

        # Test 1: Single turbine, zero covariance → should return constant mean
        CholZero = zeros(1, 1)
        wind = WindVelType(10.0, CholZero)
        vel = getWindSpeedT(Velocity_Constant_wErrorCov(), wind, [1])
        @test length(vel) == 1
        @test isapprox(vel[1], 10.0; atol=1e-6)

        # Test 2: Multiple turbines, identity covariance → should add Gaussian noise
        CholId = I(3)  # 3x3 identity matrix
        wind2 = WindVelType(12.0, Matrix{Float64}(CholId))
        vel2 = getWindSpeedT(Velocity_Constant_wErrorCov(), wind2, [1, 2, 3])
        @test length(vel2) == 3
        @test all(x -> typeof(x) <: Float64, vel2)

        # Test 3: Check statistical behavior over many samples
        samples = [getWindSpeedT(Velocity_Constant_wErrorCov(), wind, [1])[1] for _ in 1:10000]
        mean_sample = mean(samples)
        std_sample = std(samples)
        @test isapprox(mean_sample, 10.0; atol=0.1)
        @test isapprox(std_sample, 0.0; atol=0.01)  # Since CholZero ⇒ no noise

        # Test 4: Custom CholSig with known transformation
        chol_custom = [1.0 0.0; 0.5 0.5]
        wind3 = WindVelType(5.0, chol_custom)
        vel3 = getWindSpeedT(Velocity_Constant_wErrorCov(), wind3, [1, 2])
        @test length(vel3) == 2
    end
end