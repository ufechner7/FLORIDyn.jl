# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause
 
using FLORIDyn, Test, Statistics, MAT, DataFrames

matlab_file   = "test/data/input_setUpTmpWFAndRun_2_steps.mat"
after_interpolateOPs_T_file = "test/data/after_interpolateOPs_T.mat"
vars = matread(matlab_file)
vars_after_interpolateOPs_T = matread(after_interpolateOPs_T_file)
wf_dict = vars["T"]
wf_dict_ref = vars_after_interpolateOPs_T["T"]

if !isdefined(Main, :TestHelpers)
    include("test_helpers.jl")
end
using .TestHelpers

before_interpolateOPs_T_file = "test/data/before_interpolateOPs_T.mat"
vars_before_interpolateOPs_T = matread(before_interpolateOPs_T_file)
wf_dict_03 = vars_before_interpolateOPs_T["T"]

after_interpolateOPs_T_file = "test/data/after_interpolateOPs_T.mat"
vars_after_interpolateOPs_T = matread(after_interpolateOPs_T_file)
wf_dict_02 = vars_after_interpolateOPs_T["T"]

@testset "interpolateOPs_basic" begin
    global wf, wf_ref_02
    wf = convert_wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T
    wf.intOPs = interpolateOPs(wf)

    wf_ref_02 = convert_wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T
    if !compare_windFarms(wf_ref_02, wf; detailed=false, tolerance=1e-6)
        @warn "WindFarm does not match reference after interpolateOPs"
        compare_windFarms(wf_ref_02, wf; detailed=true, tolerance=1e-6)
    end
end


nothing
