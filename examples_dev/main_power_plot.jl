# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Minimal example of how to run a simulation using FLORIDyn.jl
using FLORIDyn, TerminalPager, DistributedNext, Statistics
if Threads.nthreads() == 1; using ControlPlots; end

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
include("../examples/remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions
wf = initSimulation(wf, sim)

vis.online = false
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
# plot_measurements(wf, md, vis; separated=false, msr=VelReduction, plt, pltctrl)

data_column = "ForeignReduction"
title = "Velocity Reduction"
ylabel = "Rel. Wind Speed [%]"
msr_name = "msr_velocity_reduction"

times, plot_data, turbine_labels, subplot_labels = FLORIDyn.prepare_large_plot_inputs(wf, md, data_column, ylabel; simple=true)
# plot_x(times, plot_data...; ylabels=turbine_labels, labels=subplot_labels,
#         fig=title, xlabel="rel_time [s]", ysize=9, bottom=0.02, pltctrl, legend_size=6, loc="center left")
nT = wf.nT
power_sum = zeros(length(times))
for iT in 1:nT
    rel_speed = plot_data[1][iT] ./ 100
    rel_power = rel_speed .^3
    power_sum .+= rel_power
end
power_sum ./= nT

plt.plot(times,power_sum)
plt.grid(true)
