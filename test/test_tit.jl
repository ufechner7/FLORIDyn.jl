# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test

@testset "wind turbulence" begin
    WindTi = 0.1
    iT = [1, 2, 3]
    Ti = getWindTiT(TI_Constant(), WindTi, iT)
    # Ti will be [0.1, 0.1, 0.1]
    @test length(Ti) == 3
    @test Ti[1] ≈ 0.1
    @test Ti[2] ≈ 0.1
    @test Ti[3] ≈ 0.1
    # Example WindTi data: times 0.0, 1.0, 2.0; two turbines
WindTi = [
    0.0  0.1  0.2;
    1.0  0.2  0.3;
    2.0  0.3  0.5
]

@testset "getWindTiT_EnKF Tests" begin

    # Test exact match at first time point
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 1, 0.0) ≈ 0.1
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 2, 0.0) ≈ 0.2

    # Test exact match at last time point
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 1, 2.0) ≈ 0.3
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 2, 2.0) ≈ 0.5

    # Test interpolation between points
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 1, 0.5) ≈ 0.15
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 2, 0.5) ≈ 0.25

    # Test multiple turbine indices
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, [1,2], 1.0) ≈ [0.2, 0.3]

    # Test clamping below minimum time
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 1, -1.0) ≈ 0.1
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 2, -1.0) ≈ 0.2

    # Test clamping above maximum time
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 1, 3.0) ≈ 0.3
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 2, 3.0) ≈ 0.5

    # Test with non-integer time
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 1, 1.5) ≈ 0.25
    @test getWindTiT_EnKF(TI_EnKF_InterpTurbine(), WindTi, 2, 1.5) ≈ 0.4
end

end
