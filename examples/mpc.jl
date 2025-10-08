# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Main script to run a turbine group control (TGC) simulation with FLORIDyn.jl
# using a precomputed induction matrix for feed-forward control.
# TGC shall be extended to full model predictive control (MPC) in a future example

using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext, DataFrames
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"

USE_TGC = true
USE_STEP = false
USE_FEED_FORWARD = true
ONLINE = false

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
con.yaw_fixed = 270.0
wind.input_dir="Constant"
wind.dir_fixed = 270.0
induction = calc_induction_per_group(1, 0)
set_induction!(ta, induction)

time_step = sim.time_step  # seconds
t_end = sim.end_time - sim.start_time  # relative end time in seconds
con.induction_data = calc_induction_matrix(ta, con, time_step, t_end)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
if USE_FEED_FORWARD
    set.induction_mode = Induction_TGC()
else
    set.induction_mode = Induction_Constant()
end
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
toc()


function run_simulation(set_induction::AbstractMatrix)
    global set, wind, con, floridyn, floris, sim, ta, vis 
    con.induction_data = set_induction
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    vis.online = ONLINE
    @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    # TODO calculate a power time series
    return wf, md, mi
end
induction_data = calc_induction_matrix(ta, con, time_step, t_end)
wf, md, mi = run_simulation(induction_data)

# Calculate total wind farm power by grouping by time and summing turbine powers
total_power_df = combine(groupby(md, :Time), :PowerGen => sum => :TotalPower)

# plot_rmt(md.Time, md.PowerGen)
plot_rmt(total_power_df.Time, total_power_df.TotalPower; xlabel="Time (s)", ylabel="Total Power (MW)", 
         title="Total Wind Farm Power")