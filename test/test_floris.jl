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
    @testset "getVars" begin
        # Define dummy test parameters

        RPs = [881.928  -185.932  -61.0219; 
               886.024  -170.647  -43.9671; 
               888.637  -160.896  -23.0056] # Test with 3 downstream distances
        a   = [1.0, 1.0, 1.0]                 # Dummy placeholder, not used
        C_T = [0.8844, 0.8844, 0.8844]      # Thrust coefficients
        yaw = [0.0, 0.0, 0.0]               # Yaw angles
        TI  = [0.06, 0.06, 0.06]             # Ambient intensity
        TI0 = [0.0062, 0.0062, 0.0062]      # Local turbulence
        D = 178.4                           # Rotor diameter

        # Parameter object
        struct Params
            k_a::Float64
            k_b::Float64
            alpha::Float64
            beta::Float64
        end

        param = Params(0.3837, 0.003678, 2.32, 0.154)

        sig_y, sig_z, C_T_out, x_0, delta, pc_y, pc_z = getVars(RPs, a, C_T, yaw, TI, TI0, param, D)

        @test length(sig_y) == 3
        @test length(sig_z) == 3
        @test size(delta) == (3, 2)
        @test all(x -> x ≥ 0, sig_y)
        @test all(x -> x ≥ 0, sig_z)
        @test all(x -> x ≥ 0, x_0)

        # Check for known edge values
        # @test sig_y[1] ≈ D / sqrt(8) * cos(yaw[1]) atol=1e-2
        # @test sig_z[1] ≈ D / sqrt(8) atol=1e-2

        # Check delta at RPs = 0 is only delta_nfw
        @test delta[1, 1] ≈ 0.0 atol=1e-5

        # Check pc_y and pc_z [should match D*cos(yaw), D] at RPs = 0
        # @test pc_y[1] ≈ D * cos(yaw[1]) atol=1e-4
        # @test pc_z[1] ≈ D atol=1e-4

        # Check no NaN or Inf in outputs
        for arr in (sig_y, sig_z, x_0, delta[:,1], pc_y, pc_z)
            @test all(isfinite, arr)
        end
    end
    @testset "runFLORIS" begin
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), 
                       Direction_All(), Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA())
        LocationT = [600.0  2400.0  119.0]
        States_WF = [8.2     255.0    0.062  255.0]
        States_T  = [0.33      0.0    0.06]
        D = 178.4
        paramFLORIS = FLORIDyn.Floris(
            alpha = 2.32,
            beta = 0.154,
            k_a = 0.3837,
            k_b = 0.0037,
            k_fa = 0.73,
            k_fb = 0.8325,
            k_fc = 0.0325,
            k_fd = -0.32,
            eta = 1,
            p_p = 2.2,
            airDen = 1.225,
            TIexp = 3,
            rotor_points = 50
        )
        windshear = WindShear(0.08, 1.0)
        T_red_arr, T_aTI_arr, T_Ueff, T_weight = runFLORIS(set::Settings, LocationT, States_WF, 
                                                           States_T, D, paramFLORIS, windshear)
        # Main.@infiltrate
        @test T_red_arr ≈ 0.9941836044148462
        @test isnothing(T_aTI_arr)
        @test isnothing(T_Ueff)
        @test isnothing(T_weight)
    end

end
nothing