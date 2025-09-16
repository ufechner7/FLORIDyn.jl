# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# This script is used to identify memory allocations
using FLORIDyn, DistributedNext
if Threads.nthreads() == 1; 
    @info "Running in single-threaded mode, using ControlPlots directly"
    using ControlPlots; 
end

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

vis = Vis(vis_file)
vis.show_plots = false  # Enable/disable showing plots during simulation
if (@isdefined plt) && ! isnothing(plt)
    @info "Using existing plt instance"
    plt.ion()
else
    plt = nothing
end

# Automatic parallel/threading setup
include("../examples/remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 2

vis.online = true
# Clean up any existing PNG files in video folder before starting
cleanup_video_folder()
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
nothing

# Memory allocation statistics on laptop, battery
# Baseline: 1.455825 seconds (54.85 M allocations: 5.027 GiB, 35.61% gc time)
# 219.869886 seconds (61.49 M allocations: 5.326 GiB, 0.26% gc time, 6.08% compilation time)
# Memory allocation statistics on desktop
# 360.950426 seconds (58.65 M allocations: 4.438 GiB, 0.11% gc time, 3.07% compilation time)
# 173.127742 seconds (55.33 M allocations: 4.326 GiB, 0.24% gc time, 4.82% compilation time)
# 127.343394 seconds (52.42 M allocations: 3.288 GiB, 0.30% gc time, 6.40% compilation time)
