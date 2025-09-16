# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Short testcase with online visualisation
using DistributedNext, Timers, ControlPlots, FLORIDyn

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

# Load vis settings from YAML file
vis = Vis(vis_file)
vis.online = true

# Automatic parallel/threading setup
include("remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

# prepare the simulation
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 50
wf = initSimulation(wf, sim)

cleanup_video_folder()
# Unified function automatically handles multi-threading vs single-threading
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

nothing
