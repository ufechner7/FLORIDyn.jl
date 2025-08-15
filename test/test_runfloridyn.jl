# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause
 
using FLORIDyn, Test, Statistics, MAT, DataFrames

if !isdefined(Main, :TestHelpers)
    include("test_helpers.jl")
end
using .TestHelpers

matlab_file   = "test/data/input_setUpTmpWFAndRun_2_steps.mat"
vars = matread(matlab_file)
wf_dict = vars["T"]


@testset "runfloridyn_vs_matlab" begin
    global wf1, wf_old, floris
    settings_file = "data/2021_9T_Data.yaml"
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    # create settings struct
    set = Settings(wind, sim, con)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    wf.dep = findTurbineGroups(wf, floridyn)
    wf.intOPs = interpolateOPs(wf)
    sim.n_sim_steps = 2
    wf1 = convert_wf_dict2windfarm(wf_dict)
    wf_old = deepcopy(wf)
    tmpM, wf = setUpTmpWFAndRun(set, wf1, floris, wind)
    tmpM_ref = [0, 9.094956525547340e-02, 0.1184213787906785, 0, 9.094956525547351e-02,
                0.1184213787906786, 0, 9.094956525547340e-02, 1.010889141095400e-01]
    @test all(tmpM[:,2] .≈ tmpM_ref)
    # compare_windFarms(wf1, wf_old; detailed=true)
end

nothing
