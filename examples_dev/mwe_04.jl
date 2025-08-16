# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# This script is used to identify memory allocations
using FLORIDyn, TerminalPager, DistributedNext
if Threads.nthreads() == 1; using ControlPlots; end

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
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 5

vis.online = true
# Clean up any existing PNG files in video folder before starting
cleanup_video_folder()
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
nothing

# Memory allocation statistics on laptop, battery
# Baseline: 2.654642 seconds (137.13 M allocations: 12.567 GiB, 30.04% gc time) 