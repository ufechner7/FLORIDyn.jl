# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Simplified test case plotting the flow field
using DistributedNext, Timers, ControlPlots, FLORIDyn

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

# Load vis settings from YAML file
vis = Vis(vis_file)

# Automatic setup based on threading
include("remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic threading/parallel detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 195

# Run initial conditions
wf = initSimulation(wf, sim)

vis.online = false
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
@time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)

# Single line plotting - automatically dispatches based on available resources
@time plot_flow_field(wf, X, Y, Z, vis; msr=1, plt)

nothing
