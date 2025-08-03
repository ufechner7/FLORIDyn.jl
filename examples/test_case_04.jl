# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Distributed, Timers, ControlPlots, FLORIDyn

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
PARALLEL = false
THREADING = true

if PARALLEL
    tic()
    include("../src/visualisation/remote_plotting.jl") 
    init_plotting()  # This now returns the main process plt and creates plt on workers
    toc()
end

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct
set = Settings(wind, sim, con, PARALLEL, THREADING)
# prepare the simulation
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wf = initSimulation(wf, sim)
wf, md, mi = runFLORIDyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)

@everywhere function plot_measurements(wf, md, vis; separated)
    # Create a fresh plt instance just for this task
    local_plt = ControlPlots.plt
    return plotMeasurements(local_plt, wf, md, vis; separated=separated)
end

if set.parallel
    @time @spawnat 2 plot_measurements(wf, md, vis; separated=true)
else
    plotMeasurements(plt, wf, md, vis; separated=true)
end

nothing
