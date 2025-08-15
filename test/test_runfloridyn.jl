# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause
 
using FLORIDyn, Test, Statistics, MAT, DataFrames

if !isdefined(Main, :TestHelpers)
    include("test_helpers.jl")
end
using .TestHelpers

matlab_file   = "test/data/after_prepare_simulation_T.mat"
vars = matread(matlab_file)
wf_dict = vars["T"]
matlab_file   = "test/data/after_one_step_T.mat"
vars = matread(matlab_file)
wf_dict_01 = vars["T"]

after_interpolateOPs_T_file = "test/data/after_interpolateOPs_T.mat"
vars_after_interpolateOPs_T = matread(after_interpolateOPs_T_file)
wf_dict_02 = vars_after_interpolateOPs_T["T"]

before_interpolateOPs_T_file = "test/data/before_interpolateOPs_T.mat"
vars_before_interpolateOPs_T = matread(before_interpolateOPs_T_file)
wf_dict_03 = vars_before_interpolateOPs_T["T"]


@testset "runfloridyn_vs_matlab" begin
    global wf, wf_ref, wf_ref_03, wf_debug
    settings_file = "data/2021_9T_Data.yaml"
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    # create settings struct
    set = Settings(wind, sim, con)
    vis = Vis(online=false, save=false, rel_v_min=20.0, up_int = 4)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    sim.n_sim_steps = 1
    wf_ref = wf_dict2windfarm(wf_dict) # after_prepare_simulation_T
    @test compare_windFarms(wf, wf_ref; detailed=false)
    wf_debug = [WindFarm(), WindFarm()]
    wf, md, mi = runFLORIDyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris; debug=wf_debug)

    wf_ref_03 = wf_dict2windfarm(wf_dict_03) # before_interpolateOPs_T
    if !compare_windFarms(wf_ref_03, wf_debug[2]; detailed=false, tolerance=1e-6)
        @warn "WindFarm does not match reference before interpolateOPs"
    end

    wf_ref_02 = wf_dict2windfarm(wf_dict_02) # after_interpolateOPs_T
    @test compare_windFarms(wf_ref_02, wf_debug[1]; detailed=false, tolerance=1e-6)

    wf_ref_01 = wf_dict2windfarm(wf_dict_01)
    @test compare_windFarms(wf, wf_ref_01; detailed=false, tolerance=1e-6)
    @test size(md) == (9, 6) # from Matlab
    # @test minimum(md.ForeignReduction) ≈ 72.56141032518147 # Matlab: 73.8438
    # @test mean(md.ForeignReduction)    ≈ 98.54433712619702 # Matlab: 98.
    
end

nothing
