# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Test
using FLORIDyn

# Helper to build a minimal WindFarm for tests
function _minimal_wf(nT::Int)
    # States_WF: we need at least 2 columns (velocity, direction), allocate some rows for OPs
    # Simplify: assume one OP per turbine => nOP = 1, so States_WF has nT rows
    States_WF = zeros(nT, 4) # include 4 columns so indexing States_WF[StartI,2] works
    # StartI used as matrix of indices (rangeOPs start/end). For simplicity map each turbine to its own row.
    # In existing code StartI appears indexed like wf.StartI, 2 so keep 2 columns
    StartI = hcat(collect(1:nT), collect(1:nT))
    WindFarm(
        nT = nT,
        nOP = 1,
        States_WF = States_WF,
        States_OP = zeros(nT, 4),
        States_T = zeros(nT, 3),
        posBase = zeros(2, nT),
        posNac = zeros(2, nT),
        D = fill(100.0, nT),
        StartI = StartI,
        intOPs = Vector{Matrix{Float64}}(),
        Weight = [Float64[] for _ in 1:nT],
        dep = [Int[] for _ in 1:nT],
        red_arr = zeros(nT, 1),
        Names_T = ["a", "yaw", "ti"],
        Names_WF = ["u", "dir", "ti", "phi"],
        Names_OP = ["x", "y", "z", "dw"],
    )
end

# Minimal FLORIS params (only p_p accessed in getDataVel path for I_and_I)
function _minimal_floris()
    Floris(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0, 1.0, 1.225, 0, nothing)
end

# Build minimal Settings with adjustable velocity mode
function _settings(vel_mode)
    return Settings(vel_mode, Direction_Constant(), TI_Constant(), Shear_PowerLaw(), 
                    Direction_Constant(), Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_Constant(), Induction_Constant(), false, false)
end

@testset "getDataVel branches" begin
    nT = 3
    wf = _minimal_wf(nT)
    floris = _minimal_floris()
    tmp_m = ones(nT, 1) .* 2.0 # reduction factors (so division effect visible)

    # 1. Default branch (Interpolation) -> else branch
    begin
        wind = Wind(
            "Interpolation",  # input_vel triggers else path because not matched above and vel_mode will be Velocity_Interpolation
            "Constant", "Constant", "PowerLaw",
            FLORIDyn.WindCorrection("None","None","None"),
            FLORIDyn.WindPerturbation(false,0.0,false,0.0,false,0.0),
            [0.0 8.0; 10.0 10.0], # placeholder, will be replaced below
            nothing, 0.05, nothing
        )
        # Provide simple time-speed matrix for Velocity_Interpolation
        wind.vel = [0.0 8.0; 10.0 10.0]
        set = _settings(Velocity_Interpolation())
        wf.States_WF[:,2] .= 270.0 # direction column (unused in this branch)
        u, wind2 = getDataVel(set, wind, wf, 5.0, tmp_m, floris)
        @test length(u) == nT
        @test all(u .≈ fill(9.0, nT)) # linear interpolation mid-way
        @test wind === wind2
    end

    # 2. I_and_I branch currently not fully implemented (see warning in source).
    # Provide a minimal placeholder and mark expectation as broken.
    @testset "I_and_I branch (broken)" begin
        wind = Wind(
            "I_and_I", "Constant", "Constant", "PowerLaw",
            FLORIDyn.WindCorrection("None","None","None"),
            FLORIDyn.WindPerturbation(false,0.0,false,0.0,false,0.0),
            nothing, nothing, 0.05, nothing
        )
        set = _settings(Velocity_I_and_I())
        # Expect a failure until implementation is completed
        @test_broken getDataVel(set, wind, wf, 5.0, tmp_m, floris)
    end

    # 3. RW_with_Mean branch (currently non-functional in source: call signature to getWindSpeedT has no matching method)
    # Document current behavior: expect a MethodError until implementation added.
    begin
        wind = Wind(
            "RW_with_Mean", "Constant", "Constant", "PowerLaw",
            FLORIDyn.WindCorrection("None","None","None"),
            FLORIDyn.WindPerturbation(false,0.0,false,0.0,false,0.0),
            0.0, nothing, 0.05, nothing
        )
        set = _settings(Velocity_Constant())
        wf.States_WF[:,1] .= [7.0, 8.0, 9.0]
        @test_throws MethodError getDataVel(set, wind, wf, 0.0, tmp_m, floris)
    end

    # 4. EnKF_InterpTurbine branch
    @testset "EnKF branch" begin
        wind = Wind(
            "EnKF_InterpTurbine", "Constant", "Constant", "PowerLaw",
            FLORIDyn.WindCorrection("None","None","None"),
            FLORIDyn.WindPerturbation(false,0.0,false,0.0,false,0.0),
            [0.0 6.0 7.0 8.0; 10.0 8.0 9.0 10.0], nothing, 0.05, nothing
        )
        set = _settings(Velocity_EnKF_InterpTurbine())
        # Mid-point interpolation (t = 5.0) => linear between rows
        u_mid, wind2 = getDataVel(set, wind, wf, 5.0, tmp_m, floris)
        @test wind2 === wind
        @test length(u_mid) == nT
        @test all(isapprox.(u_mid, [7.0, 8.0, 9.0]; atol=1e-8))

        # Lower bound clamp (t < first time) - just assert clamped values
        u_low, _ = getDataVel(set, wind, wf, -1.0, tmp_m, floris)
        @test all(isapprox.(u_low, [6.0, 7.0, 8.0]; atol=1e-8))

        # Upper bound clamp (t > last time) - just assert clamped values
        u_high, _ = getDataVel(set, wind, wf, 15.0, tmp_m, floris)
        @test all(isapprox.(u_high, [8.0, 9.0, 10.0]; atol=1e-8))

        # Sanity: different query times produce distinct vectors
        @test u_mid != [6.0,7.0,8.0] && u_mid != [8.0,9.0,10.0]
    end
end
