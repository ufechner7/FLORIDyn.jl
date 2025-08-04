# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Minimal example
using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext, ControlPlots
toc()

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)

# Automatic parallel/threading setup
tic()
include("../src/visualisation/smart_plotting.jl")
toc()

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions
wf = initSimulation(wf, sim)
toc()

vis.online = false
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
@time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
@time plot_flow_field(wf, X, Y, Z, vis; msr=1, plt=plt)
nothing
