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
    @testset "calcCt" begin
        # Test typical values
        @test isapprox(calcCt(0.0, nothing), 0.0)
        @test isapprox(calcCt(0.25, nothing), 4 * 0.25 * 0.75)
        @test isapprox(calcCt(0.5, nothing), 4 * 0.5 * 0.5)
        
        # Test edge cases
        @test isapprox(calcCt(1.0, nothing), 0.0)
        @test isapprox(calcCt(-0.1, nothing), 4 * -0.1 * (1 + 0.1))  # input out of typical range
        
        # Test vectorized input
        a_values = [0.0, 0.1, 0.2, 0.3]
        expected = 4 .* a_values .* (1 .- a_values)
        @test all(isapprox.(calcCt.(a_values, Ref(nothing)), expected))
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
    @testset "centerline" begin
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
    result = zeros(size(States_OP,1), 2)
    centerline!(result, States_OP, States_T, States_WF, paramFLORIS, D)

        # # Check output size
    @test size(result) == (200, 2)
    @test result == zeros(200, 2)
    end 
    @testset "getVars" begin
        # Define dummy test parameters

        RPs = [881.928  -185.932  -61.0219; 
               886.024  -170.647  -43.9671; 
               888.637  -160.896  -23.0056] # Test with 3 downstream distances
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

        floris = FLORIDyn.Floris(
            alpha = 2.32,
            beta = 0.154,
            k_a = 0.3837,
            k_b = 0.003678,
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

    # Use in-place API
    n = size(RPs, 1)
    sig_y = Vector{Float64}(undef, n)
    sig_z = Vector{Float64}(undef, n)
    x_0   = Vector{Float64}(undef, n)
    delta = Matrix{Float64}(undef, n, 2)
    pc_y  = Vector{Float64}(undef, n)
    pc_z  = Vector{Float64}(undef, n)
    getVars!(sig_y, sig_z, x_0, delta, pc_y, pc_z, RPs, C_T, yaw, TI, TI0, floris, D)
    C_T_out = C_T

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
    @testset "runFLORIS!" begin
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), 
                       Direction_All(), Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), false, false)
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
        # Create buffers for the test
        buffers = FLORIDyn.FLORISBuffers(50)  # Use 50 rotor points as defined in paramFLORIS
    runFLORIS!(buffers, set::Settings, LocationT, States_WF, 
          States_T, D, paramFLORIS, windshear)
    # Validate buffer-populated outputs
    @test buffers.T_red_arr[1] ≈ 0.9941836044148462
    @test isempty(buffers.T_aTI_arr)
    @test isempty(buffers.T_Ueff)
    @test isempty(buffers.T_weight)

        # Additional test: Check that runFLORIS! handles multiple turbines (dummy example)
        LocationT_multi = [600.0 2400.0 119.0;
                           1200.0 2600.0 119.0] 
        States_WF = [8.2  255.0  0.062  255.0;
                     8.2  255.0  0.062  255.0]
        States_T_multi = [0.33 0.0 0.06;
                          0.33 0.0 0.06]
        D = [178.4, 178.4]
        # Create buffers for multi-turbine test
        buffers_multi = FLORIDyn.FLORISBuffers(50)  # Use 50 rotor points as defined in paramFLORIS
    runFLORIS!(buffers_multi, set, LocationT_multi, States_WF, States_T_multi, D, 
          paramFLORIS, windshear)
    @test length(buffers_multi.T_red_arr) == 2
    @test length(buffers_multi.T_aTI_arr) == 1
    @test length(buffers_multi.T_Ueff) == 1
    @test length(buffers_multi.T_weight) == 1 
    end

    @testset "getUadv" begin
        # Test 1: Basic functionality with single observation point
        states_op = [881.928  -185.932  -61.0219  1000.0  0.0  0.0]  # Single OP with downstream distance 1000m
        states_t = [0.33  0.0  0.06]   # Axial induction, yaw angle (deg), turbulence intensity
        states_wf = [8.2  255.0  0.062]  # Wind speed, direction, ambient TI
        
        floris = FLORIDyn.Floris(
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
        
        d_rotor = 178.4
        
        result = getUadv(states_op, states_t, states_wf, floris, d_rotor)
        
        # Basic checks
        @test length(result) == 1
        @test 0.5 <= result[1] <= 1.0  # Advection speed should be between 50% and 100% of freestream
        @test isfinite(result[1])
        
        # Test 2: Multiple observation points with different downstream distances
        states_op_multi = [
            0.0     0.0    0.0    0.0    0.0  0.0;     # At rotor plane
            100.0   0.0    0.0    500.0  0.0  0.0;     # Before core length  
            200.0   0.0    0.0    1500.0 0.0  0.0;     # After core length
            300.0   0.0    0.0    3000.0 0.0  0.0      # Far downstream
        ]
        states_t_multi = repeat([0.33 0.0 0.06], 4, 1)
        states_wf_multi = repeat([8.2 255.0 0.062], 4, 1)
        
        result_multi = getUadv(states_op_multi, states_t_multi, states_wf_multi, floris, d_rotor)
        
        @test length(result_multi) == 4
        @test all(0.5 .<= result_multi .<= 1.0)
        @test all(isfinite.(result_multi))
        
        # Advection speed should increase with distance (wake recovers)
        @test result_multi[1] <= result_multi[2] <= result_multi[3] <= result_multi[4]
        
        # Test 3: Effect of thrust coefficient (higher CT should reduce advection speed)
        states_t_low_ct = [0.1 0.0 0.06]   # Lower axial induction -> lower CT
        states_t_high_ct = [0.4 0.0 0.06]  # Higher axial induction -> higher CT
        
        result_low_ct = getUadv(states_op, states_t_low_ct, states_wf, floris, d_rotor)
        result_high_ct = getUadv(states_op, states_t_high_ct, states_wf, floris, d_rotor)
        
        @test result_low_ct[1] > result_high_ct[1]  # Lower CT should give higher advection speed
        
        # Test 4: Effect of yaw angle
        states_t_yawed = [0.33 30.0 0.06]  # 30 degree yaw
        result_yawed = getUadv(states_op, states_t_yawed, states_wf, floris, d_rotor)
        
        @test isfinite(result_yawed[1])
        @test 0.5 <= result_yawed[1] <= 1.0
        
        # Test 5: Edge case - zero downstream distance (at rotor plane)
        states_op_zero = [0.0  0.0  0.0  0.0  0.0  0.0]
        result_zero = getUadv(states_op_zero, states_t, states_wf, floris, d_rotor)
        
        @test length(result_zero) == 1
        @test isfinite(result_zero[1])
        @test 0.5 <= result_zero[1] <= 1.0
        
        # Test 6: Effect of turbulence intensity
        states_wf_low_ti = [8.2 255.0 0.01]   # Low ambient TI
        states_wf_high_ti = [8.2 255.0 0.15]  # High ambient TI
        
        result_low_ti = getUadv(states_op, states_t, states_wf_low_ti, floris, d_rotor)
        result_high_ti = getUadv(states_op, states_t, states_wf_high_ti, floris, d_rotor)
        
        @test isfinite(result_low_ti[1])
        @test isfinite(result_high_ti[1])
        
        # Test 7: Mathematical consistency - check that formula is implemented correctly
        # For a point far downstream, we can verify the calculation manually
        states_op_far = [0.0  0.0  0.0  5000.0  0.0  0.0]  # Far downstream
        result_far = getUadv(states_op_far, states_t, states_wf, floris, d_rotor)
        
        # At far downstream, should approach the far-field formula
        # Manually calculate for comparison
        C_T_test = 4 * 0.33 * (1 - 0.33)  # calcCt(0.33, 0.0)
        yaw_test = 0.0
        I_test = sqrt(0.06^2 + 0.062^2)
        
        # Calculate core length
        x_0_test = (cos(yaw_test) * (1 + sqrt(1 - C_T_test)) / 
                   (sqrt(2) * (floris.alpha * I_test + floris.beta * (1 - sqrt(1 - C_T_test))))) * d_rotor
        
        # Since 5000m >> x_0_test, we're in far field
        k_y_test = floris.k_a * I_test + floris.k_b
        sig_y_div_D_test = max(5000.0 - x_0_test, 0.0) * k_y_test / d_rotor + min(5000.0 / x_0_test, 1.0) * cos(yaw_test) / sqrt(8)
        sig_z_div_D_test = max(5000.0 - x_0_test, 0.0) * k_y_test / d_rotor + min(5000.0 / x_0_test, 1.0) / sqrt(8)
        
        U_cen_div_U_inf_test = sqrt(max(1 - (C_T_test * cos(yaw_test)) / (8 * sig_y_div_D_test * sig_z_div_D_test), 0.0))
        expected_result = 0.5 * (1 + U_cen_div_U_inf_test)
        
        @test isapprox(result_far[1], expected_result, rtol=1e-10)
        
        # Test 8: Vector operations consistency
        # Test that broadcasting works correctly
        states_op_vector = repeat(states_op, 5, 1)
        states_t_vector = repeat(states_t, 5, 1)
        states_wf_vector = repeat(states_wf, 5, 1)
        
        result_vector = getUadv(states_op_vector, states_t_vector, states_wf_vector, floris, d_rotor)
        
        @test length(result_vector) == 5
        @test all(result_vector .≈ result[1])  # Should be identical since inputs are identical
        
        # Test 9: Input validation - check behavior with extreme values
        states_t_extreme = [0.0 0.0 0.0]  # Zero axial induction
        result_extreme = getUadv(states_op, states_t_extreme, states_wf, floris, d_rotor)
        
        @test isfinite(result_extreme[1])
        @test result_extreme[1] ≈ 1.0  # Should approach 1.0 when CT = 0
        
        # Test 10: Yaw angle effects more thoroughly
        yaw_angles = [0.0, 15.0, 30.0, 45.0]
        results_yaw = zeros(4)
        
        for (i, yaw) in enumerate(yaw_angles)
            states_t_yaw = [0.33 yaw 0.06]
            results_yaw[i] = getUadv(states_op, states_t_yaw, states_wf, floris, d_rotor)[1]
        end
        
        @test all(isfinite.(results_yaw))
        @test all(0.5 .<= results_yaw .<= 1.0)
        # Higher yaw angles should generally reduce the wake effect (increase advection speed)
        @test results_yaw[1] <= results_yaw[end]  # 0° yaw should have lower advection than 45°
    end

end
nothing