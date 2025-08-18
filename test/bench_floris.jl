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

"""
Quick correctness check: single-turbine call using preallocated buffers.
"""
# Create buffers for single-turbine case (use floris.rotor_points)
buffers_st = FLORIDyn.FLORISBuffers(paramFLORIS.rotor_points)
runFLORIS(buffers_st, set, LocationT, States_WF, 
          States_T, D, paramFLORIS, windshear)
@test buffers_st.T_red_arr[1] ≈ 0.9941836044148462
@test isempty(buffers_st.T_aTI_arr)
@test isempty(buffers_st.T_Ueff)
@test isempty(buffers_st.T_weight)

# Additional test: Check that runFLORIS handles multiple turbines (dummy example)
LocationT_multi = [600.0 2400.0 119.0;
                    1200.0 2600.0 119.0] 
States_WF = [8.2  255.0  0.062  255.0;
                8.2  255.0  0.062  255.0]
States_T_multi = [0.33 0.0 0.06;
                    0.33 0.0 0.06]
D = [178.4, 178.4]
nT = 2

println("\n=== Performance Comparison ===")

# Benchmark 1: Allocate fresh buffers each call (allocating path)
println("\n1. Allocate fresh buffers each call (allocating path):")

# Run once to get results for comparison
T_red_arr2, T_aTI_arr2, T_Ueff2, T_weight2 = begin
    # Determine rotor discretization size for buffer creation
    RPl, _ = FLORIDyn.discretizeRotor(paramFLORIS.rotor_points)
    n_points = size(RPl, 1)
    tmp_buffers = FLORIDyn.FLORISBuffers(n_points)
    runFLORIS(tmp_buffers, set, LocationT_multi, States_WF, States_T_multi, D, 
              paramFLORIS, windshear)
    tmp_buffers.T_red_arr, tmp_buffers.T_aTI_arr, tmp_buffers.T_Ueff, tmp_buffers.T_weight
end

t_wrapper = @benchmark begin
    RPl, _ = FLORIDyn.discretizeRotor($paramFLORIS.rotor_points)
    n_points = size(RPl, 1)
    tmp_buffers = FLORIDyn.FLORISBuffers(n_points)
    runFLORIS(tmp_buffers, $set, $LocationT_multi, $States_WF, $States_T_multi, $D, 
              $paramFLORIS, $windshear)
end

time_wrapper = mean(t_wrapper.times)/1e9
rel_time_wrapper = time_wrapper * 301 / 0.115  # Relative to the total time of 0.115 seconds
println("  Time: $time_wrapper seconds, relative to 0.115s: $(round(rel_time_wrapper * 100, digits=2)) %")
println("  Allocations: $(t_wrapper.allocs), Memory: $(t_wrapper.memory) bytes")

# Benchmark 2: Pre-allocated buffers (reused, zero-alloc path)
println("\n2. Pre-allocated buffers (reused, zero-alloc path):")

# Determine buffer size based on rotor discretization
if D[end] > 0
    RPl, _ = FLORIDyn.discretizeRotor(paramFLORIS.rotor_points)
    n_points = size(RPl, 1)
else
    n_points = 1
end

# Create buffers once (this is the key performance optimization)
buffers = FLORIDyn.FLORISBuffers(n_points)
println("  Created buffers for $n_points rotor discretization points")

# Run once to get results for comparison
runFLORIS(buffers, set, LocationT_multi, States_WF, States_T_multi, D, 
          paramFLORIS, windshear)
T_red_arr3, T_aTI_arr3, T_Ueff3, T_weight3 = buffers.T_red_arr, buffers.T_aTI_arr, buffers.T_Ueff, buffers.T_weight

# Benchmark with pre-allocated buffers
t_buffered = @benchmark runFLORIS(buffers, set, LocationT_multi, States_WF, States_T_multi, D, 
                                  paramFLORIS, windshear)

time_buffered = mean(t_buffered.times)/1e9
rel_time_buffered = time_buffered * 301 / 0.115
println("  Time: $time_buffered seconds, relative to 0.115s: $(round(rel_time_buffered * 100, digits=2)) %")
println("  Allocations: $(t_buffered.allocs), Memory: $(t_buffered.memory) bytes")

# Performance improvement calculation
speedup = time_wrapper / time_buffered
memory_reduction = (t_wrapper.memory - t_buffered.memory) / t_wrapper.memory * 100
alloc_reduction = (t_wrapper.allocs - t_buffered.allocs) / t_wrapper.allocs * 100

println("\n=== Performance Improvement ===")
println("  Speedup: $(round(speedup, digits=2))x faster")
println("  Memory reduction: $(round(memory_reduction, digits=1))%")
println("  Allocation reduction: $(round(alloc_reduction, digits=1))%")

# Verify results are identical
@test T_red_arr2 ≈ T_red_arr3
println("\n✓ Results are identical between allocating and preallocated versions")

println("\nRecommendation:")
println("  For single calls: Allocate buffers on the fly: FLORISBuffers(n_points) → runFLORIS(buffers, ...)")
println("  For repeated calls: Pre-allocate once and reuse: runFLORIS(buffers, set, location_t, states_wf, states_t, d_rotor, floris, windshear)")

t = t_buffered  # Use the optimized version for final reporting

time = mean(t.times)/1e9
rel_time = time * 301 / 0.115  # Relative to the total time of 0.115 seconds
println("Benchmark time: $time seconds, relative to 0.115s: $(round(rel_time * 100, digits=2)) %")
println("Allocations: $(t.allocs), in total $(t.memory) bytes.")
# Benchmark time: 1.29870531e-5 seconds, relative to 0.115s: 3.4 %
