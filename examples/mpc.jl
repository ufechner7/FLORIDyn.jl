# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Main script to run a turbine group control (TGC) simulation with FLORIDyn.jl
# using a precomputed induction matrix for feed-forward control.
# TGC shall be extended to full model predictive control (MPC) in a future example

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

# Calculate demand for each time point
time_vector = 0:time_step:t_end
demand_values = [calc_demand(t) for t in time_vector]

function run_simulation(set_induction::AbstractMatrix)
    global set, wind, con, floridyn, floris, sim, ta, vis 
    con.induction_data = set_induction
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    vis.online = ONLINE
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    # Calculate total wind farm power by grouping by time and summing turbine powers
    total_power_df = combine(groupby(md, :Time), :PowerGen => sum => :TotalPower)
end
induction_data = calc_induction_matrix(ta, con, time_step, t_end)
total_power_df = run_simulation(induction_data)

# Calculate theoretical maximum power based on turbine ratings and wind conditions
# Assumptions: free flow wind speed = 8.2 m/s, optimal axial induction factor, no yaw
wind_speed = wind.vel
a_opt = 1/3  # optimal axial induction factor (Betz limit)
Cp_opt = 4 * a_opt * (1 - a_opt)^2  # optimal power coefficient
yaw = 0.0  # no yaw angle

# Calculate maximum power for all turbines using getPower formula
# P = 0.5 * ρ * A * Cp * U^3 * η * cos(yaw)^p_p
nT = length(ta.pos[:, 1])  # number of turbines
rotor_area = π * (wf.D[1] / 2)^2  # assuming all turbines have same diameter
max_power_per_turbine = 0.5 * floris.airDen * rotor_area * Cp_opt * wind_speed^3 * floris.eta * cos(yaw)^floris.p_p / 1e6  # MW
max_power = nT * max_power_per_turbine  # total maximum power in MW

total_power_df.RelativePower = (total_power_df.TotalPower ./ max_power) .* 100  # Convert to percentage

plot_rmt(time_vector, [total_power_df.RelativePower, demand_values .* 100]; xlabel="Time [s]", xlims=(400, 1600),
         ylabel="Rel. Power Output [%]", labels=["rel_power", "rel_demand"], pltctrl)