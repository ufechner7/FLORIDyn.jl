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
if @isdefined plt
    plt.ion()
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

# With Bumper:
#  88.343513 seconds (10.41 G allocations: 825.498 GiB, 29.10% gc time, 2.40% compilation time)
# Without Bumper:
#  85.383338 seconds (11.82 G allocations: 873.189 GiB, 29.14% gc time, 1 lock conflict, 1.41% compilation time)