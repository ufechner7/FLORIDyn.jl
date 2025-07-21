# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

@testset verbose=true "setup                                      " begin
    @testset "turbineArrayProperties" begin
        tProp = turbineArrayProperties("data/2021_9T_Data.yaml")

        # Test dimensions of position matrix
        @test size(tProp.Pos)   == (9, 3)
        @test eltype(tProp.Pos) <: Real

        # Test position values of first and last turbine
        @test tProp.Pos[1, :]   == [600, 2400, 0]
        @test tProp.Pos[end, :] == [2400, 600, 0]

        # Test turbine types
        @test length(tProp.Type) == 9
        @test all(t -> t == "DTU 10MW", tProp.Type)

        # Test initial states matrix
        @test size(tProp.Init_States) == (9, 3)
        @test all(abs.(tProp.Init_States[:,1] .- 0.33) .< 1e-6)
        @test all(tProp.Init_States[:,2] .== 0)
        @test all(abs.(tProp.Init_States[:,3] .- 0.06) .< 1e-6)
    end
end