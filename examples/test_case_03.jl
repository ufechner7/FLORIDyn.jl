# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using DistributedNext, Timers, ControlPlots, FLORIDyn

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=true, save=true, rel_v_min=20.0, up_int = 4)

# Automatic parallel/threading setup
if Threads.nthreads() > 1
    tic()
    include("../src/visualisation/remote_plotting.jl") 
    init_plotting()
    toc()
end

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

# prepare the simulation
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 50
wf = initSimulation(wf, sim)

cleanup_video_folder()
# Unified function automatically handles multi-threading vs single-threading
@time wf, md, mi = smart_runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

nothing
