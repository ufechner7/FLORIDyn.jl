# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Distributed, Timers, ControlPlots, FLORIDyn

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
if Threads.nthreads() > 1
    THREADING = true
    PARALLEL  = true
else
    THREADING = false
    PARALLEL  = false
end

if PARALLEL
    tic()
    include("../src/visualisation/remote_plotting.jl") 
    init_plotting()
    toc()
end

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct
set = Settings(wind, sim, con, PARALLEL, THREADING)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 195

# Run initial conditions
wf = initSimulation(wf, sim)

vis.online = false
@time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
@time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)

# Create a completely isolated plt instance for this specific task
@everywhere function plot_flow_field(wf, X, Y, Z, vis; msr=3)
    # Create a fresh plt instance just for this task
    local_plt = ControlPlots.plt
    return plotFlowField(local_plt, wf, X, Y, Z, vis; msr=msr)
end

# Use the smart dispatcher function
@time smart_plot_flow_field(wf, X, Y, Z, vis; msr=1, plt=plt)

nothing
