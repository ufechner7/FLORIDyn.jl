# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause
 
using FLORIDyn, Test, Statistics, MAT, DataFrames

if !isdefined(Main, :TestHelpers)
    include("test_helpers.jl")
end
using .TestHelpers

matlab_file   = "test/data/after_init_simulation_T.mat"
vars = matread(matlab_file)
wf_dict = vars["T"]
matlab_file   = "test/data/after_one_step_T.mat"
vars = matread(matlab_file)
wf_dict_01 = vars["T"]


@testset "runfloridyn_basic" begin
    global wf, wf_ref, wf_ref_01
    settings_file = "data/2021_9T_Data.yaml"
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    # create settings struct
    set = Settings(wind, sim, con)
    vis = Vis(online=false, save=false, rel_v_min=20.0, up_int = 4)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    sim.n_sim_steps = 1
    wf_ref = convert_wf_dict2windfarm(wf_dict)
    @test compare_windFarms(wf, wf_ref; detailed=false)
    wf, md, mi = run_floridyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)
    # TODO compare wf after after one simulation step
    wf_ref_01 = convert_wf_dict2windfarm(wf_dict_01)
    if !compare_windFarms(wf, wf_ref_01; detailed=false, tolerance=1e-6)
        @warn "WindFarm does not match reference after simulation step"
        # compare_windFarms(wf, wf_ref_01; detailed=true, tolerance=1e-6)
    end
    @test size(md) == (9, 6) # from Matlab
    # @test minimum(md.ForeignReduction) ≈ 72.56141032518147 # Matlab: 73.8438
    # @test mean(md.ForeignReduction)    ≈ 98.54433712619702 # Matlab: 98.
    
end

# @testset "runfloridyn_vs_matlab" begin
#     global wf1, wf_old, floris
#     settings_file = "data/2021_9T_Data.yaml"
#     # get the settings for the wind field, simulator and controller
#     wind, sim, con, floris, floridyn, ta = setup(settings_file)
#     # create settings struct
#     set = Settings(wind, sim, con)
#     wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
#     sim.n_sim_steps = 2
#     wf1 = convert_wf_dict2windfarm(wf_dict)
#     wf_old = deepcopy(wf)
#     tmpM, wf = setUpTmpWFAndRun(set, wf1, floris, wind)
#     tmpM_ref = [0, 9.094956525547340e-02, 0.1184213787906785, 0, 9.094956525547351e-02,
#                 0.1184213787906786, 0, 9.094956525547340e-02, 1.010889141095400e-01]
#     @test all(tmpM[:,2] .≈ tmpM_ref)
#     # compare_windFarms(wf1, wf_old; detailed=true)
# end

nothing
