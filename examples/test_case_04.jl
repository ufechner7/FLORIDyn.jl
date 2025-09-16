# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Simplified test case plotting for the measurements
using DistributedNext, Timers, ControlPlots, FLORIDyn

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

# Load vis settings from YAML file
vis = Vis(vis_file)

# Automatic parallel/threading setup
include("remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

# prepare the simulation
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf = initSimulation(wf, sim)
wf, md, mi = run_floridyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)

# Use smart plotting function for automatic dispatch
@time plot_measurements(wf, md, vis; separated=true, plt)

nothing
