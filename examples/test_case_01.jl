# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Testcase for bug https://github.com/ufechner7/FLORIDyn.jl/issues/35
using FLORIDyn, TerminalPager, ControlPlots, DistributedNext

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

# Load vis settings from YAML file
vis = Vis(vis_file)

# Automatic parallel/threading setup
include("remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 195

# Run initial conditions
wf = initSimulation(wf, sim)

vis.online = false
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
@time Z, X, Y    = calcFlowField(set, wf, wind, floris; plt)

plot_flow_field(wf, X, Y, Z, vis; msr=EffWind, plt)
# plot_measurements(wf, md, vis; separated=true, plt)
v_min = minimum(Z[:, :, 3])
v_max = maximum(Z[:, :, 3])
println("v_min: $v_min, v_max: $v_max")
nothing
