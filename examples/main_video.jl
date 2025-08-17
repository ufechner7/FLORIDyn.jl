# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using Timers
using FLORIDyn, TerminalPager, DistributedNext
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

vis = Vis(vis_file)
vis.show_plots = true  # Enable/disable showing plots during simulation
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end

# Automatic parallel/threading setup
tic()
include("remote_plotting.jl")
toc()

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
toc()

vis.online = true
# Clean up any existing PNG files in video folder before starting
cleanup_video_folder()
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
nothing

# 61.696510 seconds (7.74 G allocations: 619.266 GiB, 29.02% gc time, 6 lock conflicts, 1.70% compilation time)
# 58.938590 seconds (7.00 G allocations: 595.492 GiB, 28.05% gc time, 7 lock conflicts, 2.15% compilation time)
# 53.673897 seconds (6.62 G allocations: 441.986 GiB, 27.28% gc time, 7 lock conflicts, 2.37% compilation time)
# 53.558160 seconds (6.56 G allocations: 440.155 GiB, 26.71% gc time, 7 lock conflicts, 2.30% compilation time)