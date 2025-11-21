# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Minimal example of how to run a simulation using FLORIDyn.jl 
# for benchmarking the 54 turbine layout.
using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext 
if Threads.nthreads() == 1; using ControlPlots; end
toc()

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"

# Load vis settings from YAML file
vis = Vis(vis_file)
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
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)
dt = 400
sim.end_time += dt
wind_dir = 270
con.yaw = "Constant"
con.yaw_data = [wind_dir;;]
wind.input_dir = "Constant"

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
toc()

vis.online = false
vis.subtitle = "54 Turbine Layout - Center-Line Model - Wind Dir: $(wind_dir)°"
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
@time Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
@time plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)
nothing
