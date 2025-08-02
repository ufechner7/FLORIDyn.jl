# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Distributed, Timers, ControlPlots, FLORIDyn

tic()
include("../src/visualisation/remote_plotting.jl") 
init_plotting()  # This now returns the main process plt and creates plt on workers
toc()

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=true, save=true, rel_v_min=20.0, up_int = 4)
PARALLEL = true
THREADING = true

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct
set = Settings(wind, sim, con, PARALLEL, THREADING)
# prepare the simulation
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf = initSimulation(wf, sim)

@everywhere function plot_flow_field(wf, X, Y, Z, vis, t_rel; msr=1)
    global plot_state
    if abs(t_rel) < 1e-6
        plot_state = nothing
    end
    local_plt = ControlPlots.plt
    plot_state = plotFlowField(plot_state, local_plt, wf, X, Y, Z, vis, t_rel; msr=msr)
    nothing
end

cleanup_video_folder()
@time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris, plot_flow_field)

# @time @spawnat 2 plot_flow_field(wf, X, Y, Z, vis; msr=3)

nothing
