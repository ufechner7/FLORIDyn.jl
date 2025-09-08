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

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"

# Load vis settings from YAML file
vis = Vis(vis_file)
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end
pltctrl = nothing
# Provide ControlPlots module only for pure sequential plotting (single-threaded, no workers)
if Threads.nthreads() == 1
    pltctrl = ControlPlots
end

# Automatic parallel/threading setup
tic()
include("remote_plotting.jl")
toc()

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions
wf = initSimulation(wf, sim)
toc()

function plot_dfs(df1, df2; fig=nothing, max_op=195)
    # Filter dataframes to keep only rows where OP <= max_op
    df1_filtered = filter(row -> row.OP <= max_op, df1)
    df2_filtered = filter(row -> row.OP <= max_op, df2)
    
    turbine_dfs1 = [df for df in groupby(df1_filtered, :Turbine)]
    turbine_dfs2 = [df for df in groupby(df2_filtered, :Turbine)]
    
    # Create vectors for each turbine plot (each containing two lines: Julia and Ref)
    turbine_plots = [[collect(turbine_dfs1[i].TI), collect(turbine_dfs2[i].TI)] for i in 1:9]
    ylabels = ["TI Turbine $i" for i in 1:9]
    labels = [["Julia", "Ref"] for _ in 1:9]

    if isnothing(fig)
        fig = "Turbine TI Comparison"
    end
    
    p = plotx(collect(turbine_dfs1[1].OP), turbine_plots...; 
              xlabel="Operating Point", ylabels=ylabels, labels=labels, ysize=10,
              fig, bottom=0.02)
    display(p)
end

vis.online = false
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
# plot_measurements(wf, md, vis; separated=false, plt)
# plotMeasurements(plt, wf, md, vis; separated=false, msr=VelReduction)    
plot_measurements(wf, md, vis; separated=false, msr=VelReduction, plt, pltctrl)
nothing
