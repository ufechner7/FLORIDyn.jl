# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Minimal example of how to run a simulation using FLORIDyn.jl 
# for benchmarking the 54 turbine layout.
using FLORIDyn

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"

# Load vis settings from YAML file
vis = Vis(vis_file)
plt = nothing

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

vis.online = false
vis.subtitle = "54 Turbine Layout - Center-Line Model - Wind Dir: $(wind_dir)°"
run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
@info "Running FLORIDyn simulation with 54 turbines"
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
nothing

# 2.2s on Ryzen 7850X
