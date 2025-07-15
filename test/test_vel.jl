# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
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
end