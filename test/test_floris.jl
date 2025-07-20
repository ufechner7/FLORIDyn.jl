# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test

@testset verbose=true "floris                                                  " begin
    @testset "discretizeRotor" begin

        # Test 1: Check output sizes
        nRP = 27  # example input
        RPs, w = discretizeRotor(nRP)
        N1 = 3
        n = round(Int, sqrt(nRP / N1))
        expected_nC = N1 * n^2

        @test size(RPs) == (expected_nC, 3)
        @test length(w) == expected_nC

        # Test 2: Check that RPs values in the 2nd and 3rd columns are in [-0.5, 0.5]
        @test all(abs.(RPs[:,2]) .<= 0.5)
        @test all(abs.(RPs[:,3]) .<= 0.5)

        # Test 3: Check that first column of RPs is all zeros (matches MATLAB code)
        @test all(RPs[:,1] .== 0.0)

        # Test 4: Check that weights sum to approximately 1
        @test isapprox(sum(w), 1.0; atol=1e-14)

        # Test 5: For different input, check outputs are plausible
        nRP2 = 48
        RPs2, w2 = discretizeRotor(nRP2)
        n2 = round(Int, sqrt(nRP2 / N1))
        @test size(RPs2) == (N1*n2^2, 3)
        @test all(abs.(RPs2[:,2]) .<= 0.5)
        @test all(abs.(RPs2[:,3]) .<= 0.5)
        @test isapprox(sum(w2), 1.0; atol=1e-14)
    end
    @testset "CalcCt" begin
        # Test typical values
        @test isapprox(CalcCt(0.0, nothing), 0.0)
        @test isapprox(CalcCt(0.25, nothing), 4 * 0.25 * 0.75)
        @test isapprox(CalcCt(0.5, nothing), 4 * 0.5 * 0.5)
        
        # Test edge cases
        @test isapprox(CalcCt(1.0, nothing), 0.0)
        @test isapprox(CalcCt(-0.1, nothing), 4 * -0.1 * (1 + 0.1))  # input out of typical range
        
        # Test vectorized input
        a_values = [0.0, 0.1, 0.2, 0.3]
        expected = 4 .* a_values .* (1 .- a_values)
        @test all(isapprox.(CalcCt.(a_values, Ref(nothing)), expected))
    end
    @testset "States() constructor" begin
        s = States()

        # Test type
        @test isa(s, States)

        # Test Turbine states
        @test s.T_names == ["a", "yaw", "TI"]
        @test s.Turbine == 3

        # Test Observation Point states
        @test s.OP_names == ["x0", "y0", "z0", "x1", "y1", "z1"]
        @test s.OP == 6

        # Test Wind Field states
        @test s.WF_names == ["wind_vel", "wind_dir", "TI0"]
        @test s.WF == 3
    end
    @testset "Centerline function" begin
        # Define dummy parameters
        States_OP = zeros(200, 6)
        States_OP[:, 4] = (0.0:0.0328:6.5272) * 1e3
        States_T = zeros(200, 3)
        States_T[:, 1] .= 0.33
        States_T[:, 3] .= 0.06 # Example turbulence intensity
        States_WF = zeros(200, 4)
        States_WF[:, 4] .= 255
        States_WF[:, 3] .= 0.0620 
        States_WF[:, 2] .= 255
        States_WF[:, 1] .= 8.2

        paramFLORIS = FLORIDyn.Floris(
            2.32,      # alpha
            0.154,     # beta
            0.3837,    # k_a
            0.0037,    # k_b
            0.73,      # k_fa
            0.8325,    # k_fb
            0.0325,    # k_fc
            -0.32,     # k_fd
            1,         # eta
            2.2,       # p_p
            1.225,     # airDen
            3,         # TIexp
            nothing
        )

        D = 178.4000 
        result = Centerline(States_OP, States_T, States_WF, paramFLORIS, D)

        # # Check output size
        @test size(result) == (200, 2)
        @test result == zeros(200, 2)
    end
end
nothing