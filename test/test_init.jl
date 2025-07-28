# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

@testset verbose=true "setup                                      " begin
    @testset "turbineArrayProperties" begin
        tProp = turbineArrayProperties("data/2021_9T_Data.yaml")

        # Test dimensions of position matrix
        @test size(tProp.pos)   == (9, 3)
        @test eltype(tProp.pos) <: Real

        # Test position values of first and last turbine
        @test tProp.pos[1, :]   == [600, 2400, 0]
        @test tProp.pos[end, :] == [2400, 600, 0]

        # Test turbine types
        @test length(tProp.type) == 9
        @test all(t -> t == "DTU 10MW", tProp.type)

        # Test initial states matrix
        @test size(tProp.init_States) == (9, 3)
        @test all(abs.(tProp.init_States[:,1] .- 0.33) .< 1e-6)
        @test all(tProp.init_States[:,2] .== 0)
        @test all(abs.(tProp.init_States[:,3] .- 0.06) .< 1e-6)

        # Test return type
        @test isa(tProp, TurbineArray)
        
        # Test field types
        @test isa(tProp.pos, Matrix{Float64})
        @test isa(tProp.type, Vector{String})
        @test isa(tProp.init_States, Matrix{Float64})
        
        # Test consistency between fields
        @test size(tProp.pos, 1) == length(tProp.type) == size(tProp.init_States, 1)
        @test size(tProp.pos, 2) == 3  # x, y, z coordinates
        @test size(tProp.init_States, 2) == 3  # a, yaw, ti
    end
end