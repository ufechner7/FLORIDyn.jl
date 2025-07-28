# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test

@testset verbose=true "wind turbulence                                         " begin
     @testset "getWindTiT(TI_Constant(), ...)" begin
        turb_mode = TI_Constant()
        wind_ti = 0.1
        iT = [1, 2, 3]
        Ti = getWindTiT(turb_mode, wind_ti, iT, nothing)
        # Ti will be [0.1, 0.1, 0.1]
        @test length(Ti) == 3
        @test Ti[1] ≈ 0.1
        @test Ti[2] ≈ 0.1
        @test Ti[3] ≈ 0.1
    end
    # Example wind_ti data: times 0.0, 1.0, 2.0; two turbines
    wind_ti = [
        0.0  0.1  0.2;
        1.0  0.2  0.3;
        2.0  0.3  0.5
    ]

        @testset "getWindTiT(TI_InterpTurbine(), ...)" begin
        # Example wind TI data: columns are [time, TI_T1, TI_T2]
        wind_ti = [
            0.0  0.10  0.20;
            1.0  0.15  0.25;
            2.0  0.20  0.30
        ]

        turb_mode = TI_InterpTurbine()

        # Test interpolation at t=0.5 for turbine 1
        @test getWindTiT(turb_mode, wind_ti, 1, 0.5) ≈ 0.125 atol=1e-8
        # Test interpolation at t=1.5 for turbine 2
        @test getWindTiT(turb_mode, wind_ti, 2, 1.5) ≈ 0.275 atol=1e-8

        # Test clamping below range: t = -1.0
        @test getWindTiT(turb_mode, wind_ti, 1, -1.0) ≈ 0.10 atol=1e-8
        # Test clamping above range: t = 3.0
        @test getWindTiT(turb_mode, wind_ti, 2, 3.0) ≈ 0.30 atol=1e-8

        # Test exact time point
        @test getWindTiT(turb_mode, wind_ti, 2, 1.0) ≈ 0.25 atol=1e-8
    end

        @testset "getWindTiT(TI_Interpolation(), ...)" begin
        dir_mode = TI_Interpolation()

        wind_ti = [
            0.0  0.10;
            10.0 0.20;
            20.0 0.30
        ]

        # Test 1: Interpolation within range
        iT = [1, 2, 3]
        t = 5.0
        Ti = getWindTiT(dir_mode, wind_ti, iT, t)
        @test Ti ≈ fill(0.15, 3)   # (0.10 + (0.20-0.10)*5/10) = 0.15

        # Test 2: Exact time match
        t2 = 10.0
        Ti2 = getWindTiT(dir_mode, wind_ti, iT, t2)
        @test Ti2 ≈ fill(0.20, 3)

        # Test 3: Time before first entry (should clamp and warn)
        t3 = -5.0
        Ti3 = getWindTiT(dir_mode, wind_ti, iT, t3)
        @test Ti3 ≈ fill(0.10, 3)

        # Test 4: Time after last entry (should clamp and warn)
        t4 = 25.0
        Ti4 = getWindTiT(dir_mode, wind_ti, iT, t4)
        @test Ti4 ≈ fill(0.30, 3)

        # Test 5: Single turbine index (scalar)
        iT5 = [2]
        t5 = 15.0
        Ti5 = getWindTiT(dir_mode, wind_ti, iT5, t5)
        @test Ti5 ≈ fill(0.25, 1)  # (0.20 + (0.30-0.20)*5/10) = 0.25
    end
end
