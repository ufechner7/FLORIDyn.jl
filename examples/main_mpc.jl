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
include("remote_plotting.jl")
include("calc_induction_matrix.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)
sim.end_time += 420
con.yaw="Constant"
wind.input_dir="Constant"

time_step = sim.time_step  # seconds
t_end = sim.end_time - sim.start_time  # relative end time in seconds
con.induction_data = calc_induction_matrix(ta, con, time_step, t_end)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
wind.dir=[270.0;;]
toc()

vis.online = false
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
# @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
# @time plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)

data_column = "ForeignReduction"
ylabel = "Rel. Wind Speed [%]"

times, plot_data, turbine_labels, subplot_labels = FLORIDyn.prepare_large_plot_inputs(wf, md, data_column, ylabel; simple=true)
nT = wf.nT
rel_power = zeros(length(times))
for iT in 1:nT
    rel_speed = plot_data[1][iT] ./ 100
    rel_power .+= rel_speed .^3
end
rel_power ./= nT
plot_rmt(times, rel_power .* 100; xlabel="Time [s]", ylabel="Rel. Power Output [%]", pltctrl)
