# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, BenchmarkTools, Parameters


set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), 
                Direction_All(), Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(),
                false, false)
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
# allocate placeholder for benchmarking allocations (not used here)
alloc = nothing
# single turbine case: default alloc, returns nothing for TI, Ueff, weight
T_red_arr, T_aTI_arr, T_Ueff, T_weight = runFLORIS(set, LocationT, States_WF, 
                                                    States_T, D, paramFLORIS, windshear)
@test T_red_arr ≈ 0.9941836044148462
@test isnothing(T_aTI_arr)
@test isnothing(T_Ueff)
@test isnothing(T_weight)

# Additional test: Check that runFLORIS handles multiple turbines (dummy example)
LocationT_multi = [600.0 2400.0 119.0;
                    1200.0 2600.0 119.0] 
nStates_WF = [8.2  255.0  0.062  255.0;
               8.2  255.0  0.062  255.0]
nStates_T = [0.33 0.0 0.06;
             0.33 0.0 0.06]
D = [178.4, 178.4]
nT = 2

T_red_arr, T_aTI_arr, T_Ueff, T_weight = runFLORIS(set, LocationT_multi, nStates_WF, 
                                                   nStates_T, D, paramFLORIS, windshear)
println(T_aTI_arr)
@test sum(T_aTI_arr) > 1e-6
nothing