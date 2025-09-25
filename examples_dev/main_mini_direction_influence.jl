# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Minimal example of how to run a simulation using FLORIDyn.jl
using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext 
if Threads.nthreads() == 1; using ControlPlots; end
toc()

settings_file, vis_file = get_default_project()[2:3]

# Load vis settings from YAML file
vis = Vis(vis_file)
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end

# Automatic parallel/threading setup
tic()
include("../examples/remote_plotting.jl")
toc()

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
set.cor_dir_mode=Direction_Influence()

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions
wf = initSimulation(wf, sim)
toc()

vis.online = false
vis.v_min = 2.0
vis.v_max = 10.0
# Matlab: 0.43s for 9T on desktop
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
@time Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
@time plot_flow_field(wf, X, Y, Z, vis; msr=EffWind, plt)
nothing
