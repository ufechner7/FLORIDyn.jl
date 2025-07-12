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

@testset "getWindTiT" begin
    # Example wind TI data: columns are [time, TI_T1, TI_T2]
    WindTi = [
        0.0  0.10  0.20;
        1.0  0.15  0.25;
        2.0  0.20  0.30
    ]

    ti_interp = TI_InterpTurbine()

    # Test interpolation at t=0.5 for turbine 1
    @test getWindTiT(ti_interp, WindTi, 1, 0.5) ≈ 0.125 atol=1e-8
    # Test interpolation at t=1.5 for turbine 2
    @test getWindTiT(ti_interp, WindTi, 2, 1.5) ≈ 0.275 atol=1e-8

    # Test clamping below range: t = -1.0
    @test getWindTiT(ti_interp, WindTi, 1, -1.0) ≈ 0.10 atol=1e-8
    # Test clamping above range: t = 3.0
    @test getWindTiT(ti_interp, WindTi, 2, 3.0) ≈ 0.30 atol=1e-8

    # Test exact time point
    @test getWindTiT(ti_interp, WindTi, 2, 1.0) ≈ 0.25 atol=1e-8
end

end
