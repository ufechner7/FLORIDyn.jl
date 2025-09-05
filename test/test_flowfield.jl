# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Testcase for bug https://github.com/ufechner7/FLORIDyn.jl/issues/35
using FLORIDyn, MAT, Test, DistributedNext, LinearAlgebra, Statistics, DataFrames

if !isdefined(Main, :TestHelpers)
    include("test_helpers.jl")
end
using .TestHelpers

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"
input_file    = "test/data/input_T_195_steps.mat"
T_ref = matread(input_file)["T"]
turbines_ref = (turbines(T_ref))  # Creates a 1800×3 DataFrame with turbine states
matlab_file   = "test/data/flowfield_xyz_195_steps.mat"
vars = matread(matlab_file)
X_ref = vars["X"]
Y_ref = vars["Y"]
Z_ref = vars["Z"]

function rel_err(a, b)
    return norm(a - b) / norm(b)
end

# Load vis settings from YAML file
vis = Vis(vis_file)

# Automatic parallel/threading setup
if !isdefined(Main, :init_plotting)
    include("../examples/remote_plotting.jl")
end

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

# create settings struct
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 195

# Run initial conditions
wf = initSimulation(wf, sim)

vis.online = false
wf, md, mi = run_floridyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)

turbines_wf = wf.turbines

@testset verbose=true "Flow Field Comparison" begin

    df1, df2 = compare_dataframes(turbines_wf, turbines_ref)
    @test size(df1, 1) == 0

    Z, X, Y    = calcFlowField(set, wf, wind, floris)
    msr = 1
    A = Z_ref[:,:,msr]
    B = Z[:,:,msr]
    A = Z_ref[:,:,msr]
    @test rel_err(A, B) < 0.001
end

nothing
