# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, BenchmarkTools, Parameters

@with_kw mutable struct AllocsF
    n::Int64 = 0      # number of floris calls
    expr1::Int64 = 0  # allocated memory of first expression
    expr2::Int64 = 0  # allocated memory of first expression
    expr3::Int64 = 0  # allocated memory of first expression
    for1::Int64 = 0   # allocated memory outer for loop
    for2::Int64 = 0   # allocated memory first inner for loop
    for3::Int64 = 0   # allocated memory first inner for loop
    if1::Int64 = 0    # allocated memory first if clause
    if2::Int64 = 0    # allocated memory first if clause
    if3::Int64 = 0    # allocated memory first if clause
end

function Base.show(io::IO, allocs::AllocsF)
    println(io, "Allocations:")
    for field_name in fieldnames(typeof(allocs))
        value = getfield(allocs, field_name)
        n = getfield(allocs, :n)
        if value > 5 && ! (field_name in [:n, :m])
            kb_value = value / 1024/ n
            println(io, "  $field_name: $(round(kb_value, digits=3)) KiB")
        end
    end
end

alloc = AllocsF()

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
T_red_arr, T_aTI_arr, T_Ueff, T_weight = runFLORIS(set::Settings, LocationT, States_WF, 
                                                    States_T, D, paramFLORIS, windshear; alloc)
# Main.@infiltrate
@test T_red_arr ≈ 0.9941836044148462
@test isnothing(T_aTI_arr)
@test isnothing(T_Ueff)
@test isnothing(T_weight)

# Additional test: Check that runFLORIS handles multiple turbines (dummy example)
LocationT_multi = [600.0 2400.0 119.0;
                    1200.0 2600.0 119.0] 
States_WF = [8.2  255.0  0.062  255.0;
                8.2  255.0  0.062  255.0]
States_T_multi = [0.33 0.0 0.06;
                    0.33 0.0 0.06]
D = [178.4, 178.4]
nT = 2

t = @benchmark T_red_arr2, T_aTI_arr2, T_Ueff2, T_weight2 = runFLORIS(set, LocationT_multi, States_WF, States_T_multi, D, 
                                                        paramFLORIS, windshear; alloc)

time = mean(t.times)/1e9
rel_time = time * 301 / 0.115  # Relative to the total time of 0.115 seconds
println("Benchmark time: $time seconds, relative to 0.115s: $(round(rel_time * 100, digits=2)) %")
println("Allocations: $(t.allocs), in total $(t.memory) bytes.")
# Benchmark time: 1.29870531e-5 seconds, relative to 0.115s: 3.4 %

println(alloc)