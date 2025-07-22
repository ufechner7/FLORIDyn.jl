# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using Random

@testset verbose=true "wind shear                                              " begin
    @testset "getWindShearT(Shear_Interpolation(), ...)" begin
        shear_mode = Shear_Interpolation()
        wind_shear = [10 0.8; 20 0.9; 30 1.0]  # Example data
        z = [5, 15, 25, 35]                   # Heights to interpolate at
        shear = getWindShearT(shear_mode, wind_shear, z)
        @test shear[1] ≈ 0.8
        @test shear[2] ≈ 0.85
        @test shear[3] ≈ 0.95
        @test shear[4] ≈ 1.0
    end
    @testset "getWindShearT(Shear_PowerLaw(), ...)" begin
        shear_mode = Shear_PowerLaw()
        # Test with scalar z_norm
        ws = WindShear(0.14, 10.0)
        z = 2.0
        expected = z^ws.alpha
        @test getWindShearT(shear_mode, ws, z) ≈ expected

        # Test with array z_norm
        z_arr = [1.0, 2.0, 4.0]
        expected_arr = z_arr .^ ws.alpha
        @test all(getWindShearT(shear_mode, ws, z_arr) .≈ expected_arr)

        # Test with alpha = 0 (should return ones)
        ws_zero = WindShear(0.0, 10.0)
        @test all(getWindShearT(shear_mode, ws_zero, z_arr) .≈ ones(length(z_arr)))

        # Test with z_norm = 0 (should return 0 unless alpha==0)
        @test getWindShearT(shear_mode, ws, 0.0) == 0.0
    end
end

nothing