# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test
@testset verbose=true "floridyncl" begin
    include("test_prepare_simulation.jl")
    @testset "angSOWFA2world" begin
        # Test 1: deg_SOWFA = 0 should give rad_World = deg2rad(270)
        @test isapprox(angSOWFA2world(0), deg2rad(270))

        # Test 2: deg_SOWFA = 90 should give rad_World = deg2rad(180)
        @test isapprox(angSOWFA2world(90), deg2rad(180))

        # Test 3: deg_SOWFA = 180 should give rad_World = deg2rad(90)
        @test isapprox(angSOWFA2world(180), deg2rad(90))

        # Test 4: deg_SOWFA = 270 should give rad_World = deg2rad(0)
        @test isapprox(angSOWFA2world(270), deg2rad(0))

        # Test 5: deg_SOWFA = 360 should give rad_World = deg2rad(-90)
        @test isapprox(angSOWFA2world(360), deg2rad(-90))
    end
end
nothing