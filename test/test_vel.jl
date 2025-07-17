# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using Statistics
using Random
using Logging

@testset verbose=true "wind velocity  " begin
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
    @testset "getWindSpeedT(Velocity_Constant_wErrorCov(), ...)" begin
        # Set fixed random seed for consistent results
        # Random.seed!(1234)

        vel_mode = Velocity_Constant_wErrorCov()

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
    WindVel = [
        0.0   6.0   8.0;
        5.0   10.0  12.0;
        10.0  14.0  18.0
    ]

    # Reference to the function for clarity
    getWS = getWindSpeedT_EnKF

    @testset "getWindSpeedT_EnKF(Velocity_EnKF_InterpTurbine(), ...)" begin

        let vel_mode = Velocity_EnKF_InterpTurbine()

            @testset "Exact time matches" begin
                @test getWS(vel_mode, WindVel, 1, 0.0) ≈ 6.0
                @test getWS(vel_mode, WindVel, 2, 10.0) ≈ 18.0
            end

            @testset "Interpolated values at midpoints" begin
                @test getWS(vel_mode, WindVel, 1, 2.5) ≈ 8.0  # halfway between 6.0 and 10.0
                @test getWS(vel_mode, WindVel, 2, 7.5) ≈ 15.0 # halfway between 12.0 and 18.0
            end

            @testset "Multiple turbine indices" begin
                result = getWS(vel_mode, WindVel, [1, 2], 2.5)
                @test result ≈ [8.0, 10.0]
            end

            @testset "Out-of-bounds: before first timestamp" begin
                buff = IOBuffer()
                with_logger(ConsoleLogger(buff)) do
                    result = getWS(vel_mode, WindVel, 1, -3.0)
                    @test result ≈ 6.0  # clamped to t = 0.0
                end
                seekstart(buff)
                log_output = String(take!(buff))
                @test occursin("out of bounds", log_output)
            end

            @testset "Out-of-bounds: after last timestamp" begin
                logbuffer = IOBuffer()
                with_logger(ConsoleLogger(logbuffer)) do
                    result = getWS(vel_mode, WindVel, 2, 15.0)
                    @test result ≈ 18.0  # clamped to t = 10.0
                end
                seekstart(logbuffer)
                logs = String(take!(logbuffer))
                @test occursin("out of bounds", logs)
            end

        end
    end
    @testset "getWindSpeedT(Velocity_InterpTurbine(), ...)" begin
        vel_mode = Velocity_InterpTurbine()
        
        # Sample data: each row is [time, turbine_1_speed, turbine_2_speed]
        WindVel = [
            0.0  5.0  6.0;
            1.0  7.0  8.0;
            2.0  9.0 10.0
        ]

        @testset "Midpoint interpolation" begin
            # At time 0.5 → linear interpolation between times 0.0 and 1.0
            result = getWindSpeedT(vel_mode, WindVel, 1, 0.5)
            @test isapprox(result, 6.0, atol=1e-8)
            
            result2 = getWindSpeedT(vel_mode, WindVel, 2, 0.5)
            @test isapprox(result2, 7.0, atol=1e-8)
        end

        @testset "Exact match on time points" begin
            @test getWindSpeedT(vel_mode, WindVel, 1, 0.0) == 5.0
            @test getWindSpeedT(vel_mode, WindVel, 2, 1.0) == 8.0
            @test getWindSpeedT(vel_mode, WindVel, 1, 2.0) == 9.0
        end

        @testset "Time before start (clamp lower bound)" begin
            result = getWindSpeedT(vel_mode, WindVel, 1, -1.0)
            @test result == 5.0
        end

        @testset "Time after end (clamp upper bound)" begin
            result = getWindSpeedT(vel_mode, WindVel, 2, 5.0)
            @test result == 10.0
        end

        @testset "Multiple index types" begin
            t = 1.5
            expected = [8.0, 9.0]
            for (i, expect) in enumerate(expected)
                @test isapprox(getWindSpeedT(vel_mode, WindVel, i, t), expect, atol=1e-8)
            end
        end
    end

    @testset "getWindSpeedT(Velocity_Interpolation(), ...)" begin
        # Example wind data (time, wind_speed)
        WindVel = [0.0 5.0; 10.0 10.0; 20.0 15.0]
        interp = Velocity_Interpolation()

        # Test 1: Exact time match
        @test getWindSpeedT(interp, WindVel, [1], 10.0) == [10.0]

        # Test 2: Interpolation between 10.0 and 20.0 (expected 12.5 => speed = 11.25)
        result = getWindSpeedT(interp, WindVel, [1,2,3], 12.5)
        @test all(x -> isapprox(x, 11.25, atol=1e-6), result)

        # Test 3: Lower bound clamp
        result_low = getWindSpeedT(interp, WindVel, [1,2], -5.0)
        @test result_low == [5.0, 5.0]

        # Test 4: Upper bound clamp
        result_high = getWindSpeedT(interp, WindVel, [1], 30.0)
        @test result_high == [15.0]

        # Test 5: Empty turbine index list
        result_empty = getWindSpeedT(interp, WindVel, Int[], 10.0)
        @test result_empty == []

        # Test 6: Single turbine index
        result_single = getWindSpeedT(interp, WindVel, [42], 15.0)
        @test result_single == [12.5]
    end

    @testset "getWindSpeedT(Velocity_Interpolation_wErrorCov(), ...)" begin
        Random.seed!(42)  # Ensure reproducibility
        
        # Simple linear data: speed = 5 at t=0, speed = 15 at t=10
        times = [0.0, 10.0]
        speeds = [5.0, 15.0]
        data = hcat(times, speeds)

        # Simple 2x2 identity Cholesky: noise still generated but uncorrelated
        chol = cholesky(Matrix(1.0I, 2, 2)).L

        wind = WindVelMatrix(data, chol)
        method = Velocity_Interpolation_wErrorCov()

        # Test 1: Midpoint interpolation
        iT = [1, 2]
        t = 5.0
        vel = getWindSpeedT(method, wind, iT, t)
        @test length(vel) == length(iT)
        @test isapprox(mean(vel), 10.0; atol=3.0)  # Within noise bounds

        # Test 2: Out-of-bounds time (before)
        t2 = -5.0
        vel2 = getWindSpeedT(method, wind, iT, t2)
        @test isapprox(mean(vel2), 5.0; atol=3.0)

        # Test 3: Out-of-bounds time (after)
        t3 = 20.0
        vel3 = getWindSpeedT(method, wind, iT, t3)
        @test isapprox(mean(vel3), 15.0; atol=3.0)

        # Test 4: Check that returned values differ (confirming noise)
        vels = [getWindSpeedT(method, wind, iT, 5.0) for _ in 1:10]
        means = [mean(v) for v in vels]
        @test std(means) > 0.0  # Some variability due to random noise

        # # Test 5: Single turbine (scalar iT input)
        # iT_scalar = 1
        # vel_scalar = getWindSpeedT(method, wind, iT_scalar, 5.0)
        # @test length(vel_scalar) == 1
        # @test isapprox(vel_scalar[1], 10.0; atol=3.0)
    end
    @testset "getWindSpeedT(Velocity_InterpTurbine_wErrorCov(), ...)" begin
        vel_mode = Velocity_InterpTurbine_wErrorCov()
        
        # Fixed Random Seed
        RNG = MersenneTwister(1234)

        # Define mock data: times = [0.0, 10.0, 20.0], wind speed at 2 turbines
        Data = [
            0.0   5.0  6.0;
            10.0 10.0 12.0;
            20.0 15.0 18.0
        ]
        
        # Cholesky factor: Identity matrix for no coupling, deterministic noise
        CholSig_zero = zeros(2,2)
        CholSig_identity = Matrix{Float64}(I, 2, 2)

        #####################################################################
        @testset "Interpolation without noise" begin
            WindVel = WindVelMatrix(Data, CholSig_zero)
            t_test = 5.0
            iT = [1, 2]

            vel = getWindSpeedT(vel_mode, WindVel, iT, t_test)
            expected = [7.5, 9.0]  # Linear interpolation at t=5.0
            @test isapprox(vel, expected, atol=1e-10)
        end

        #####################################################################
        # @testset "Clamping t below and above range" begin
        #     WindVel = WindVelMatrix(Data, CholSig_zero)

        #     # Below range
        #     vel_lo = getWindSpeedT(vel_mode, WindVel, [1], -5.0)
        #     @test isapprox(vel_lo[1], 5.0)

        #     # Above range
        #     vel_hi = getWindSpeedT(vel_mode, WindVel, [2], 30.0)
        #     @test isapprox(vel_hi[1], 18.0)
        # end

        #####################################################################
        # @testset "Noise is applied when CholSig is identity" begin
        #     WindVel = WindVelMatrix(Data, CholSig_identity)

        #     t_test = 10.0
        #     iT = [1, 2]

        #     # Use fixed seed so we can predict noise
        #     Random.seed!(RNG)
        #     vel1 = getWindSpeedT(vel_mode, WindVel, iT, t_test)

        #     Random.seed!(RNG)  # Reset RNG to get same result
        #     vel2 = getWindSpeedT(vel_mode, WindVel, iT, t_test)

        #     @test vel1 == vel2        # Result is repeatable with fixed RNG
        #     @test vel1 != [10.0, 12.0]  # Should differ from ground truth due to added noise
        # end

        #####################################################################
        # @testset "Single index access works" begin
        #     WindVel = WindVelMatrix(Data, CholSig_zero)
        #     t_test = 10.0
        #     vel = getWindSpeedT(vel_mode, WindVel, 1, t_test)
        #     @test isapprox(vel[1], 10.0)
        # end

        #####################################################################
        @testset "Interpolation at exact time point" begin
            WindVel = WindVelMatrix(Data, CholSig_zero)

            vel = getWindSpeedT(vel_mode, WindVel, [1,2], 0.0)
            @test isequal(vel, [5.0, 6.0])

            vel2 = getWindSpeedT(vel_mode, WindVel, [1,2], 20.0)
            @test isequal(vel2, [15.0, 18.0])
        end
        
    end
    @testset "getWindSpeedT(Velocity_ZOH_wErrorCov(), ...)" begin
        @testset "getWindSpeedT basic behavior" begin
            # Setup: make inputs deterministic for testing
            Vel = [1.0, 2.0]
            WindVelCholSig = Matrix{Float64}(I, 2, 2)  # identity => noise stays uncorrelated
            
            # Copy initial
            Vel_copy = copy(Vel)
            
            # Test with deterministic random seed
            FLORIDyn.set_rng(MersenneTwister(1234))
            noise = randn(FLORIDyn.RNG, 2)
            expected = Vel_copy .+ WindVelCholSig * noise
            
            # Call actual function with the same seed
            FLORIDyn.set_rng(MersenneTwister(1234))
            getWindSpeedT(Velocity_ZOH_wErrorCov(), Vel, WindVelCholSig)
            @test isapprox(Vel, expected, atol=1e-10)
        end

        @testset "getWindSpeedT zero noise" begin
            Vel = [3.0, -1.0, 7.0]
            WindVelCholSig = zeros(3, 3)  # zero matrix, noise always zero
            Vel_copy = copy(Vel)
            getWindSpeedT(Velocity_ZOH_wErrorCov(), Vel, WindVelCholSig)
            @test Vel == Vel_copy
        end

        @testset "getWindSpeedT size mismatch throws" begin
            Vel = [1.0, 2.0]
            WindVelCholSig = Matrix{Float64}(I, 3, 3)
            @test_throws DimensionMismatch getWindSpeedT(Velocity_ZOH_wErrorCov(), Vel, WindVelCholSig)
        end
    end
end