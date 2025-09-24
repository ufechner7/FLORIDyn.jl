# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Main script to run a model predictive control (MPC) simulation with FLORIDyn.jl
# using a precomputed induction matrix for feed-forward control.
using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext 
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"

USE_MPC = true
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
    set.induction_mode = Induction_MPC()
else
    set.induction_mode = Induction_Constant()
end
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
toc()

vis.online = ONLINE
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
# @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
# @time plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)

data_column = "ForeignReduction"
ylabel      = "Rel. Wind Speed [%]"

times, plot_data, turbine_labels, subplot_labels = FLORIDyn.prepare_large_plot_inputs(wf, md, data_column, ylabel; simple=true)
nT = wf.nT
rel_power = zeros(length(times))

induction_factors = zeros(nT, length(times))
for (i, sim_time) in pairs(times)
    induction_factors[:, i] = getInduction(set.induction_mode, con, (1:nT), sim_time)
end
for iT in 1:nT
    rel_speed = plot_data[1][iT] ./ 100
    induction_vec = induction_factors[iT, :]
    cp_vec = 4 * induction_vec .* (1 .- induction_vec).^2
    rel_power .+= rel_speed .^3 .* cp_vec ./ cp_max
end
rel_power ./= nT

time_vector = 0:time_step:t_end

# Calculate demand for each time point
demand_values = [calc_demand(t) for t in time_vector]

plot_rmt(times, [rel_power .* 100, demand_values .* 100]; xlabel="Time [s]", xlims=(400, 1600),
         ylabel="Rel. Power Output [%]", labels=["rel_power", "rel_demand"], pltctrl)

# Calculate Mean Square Error between rel_power and demand_values
mse = sum((rel_power[101:end] .- demand_values[101:end]).^2) / length(rel_power[101:end])
println("Root Mean Square Error (RMSE): $(round(sqrt(mse) * 100, digits=2))%")
println("Max Absolute Error:            $(round(maximum(abs.(rel_power[101:end] .- demand_values[101:end])) * 100, digits=2))%")

# Print wind conditions
println("\n--- Wind Conditions ---")
if hasfield(typeof(wind), :vel) && !isnothing(wind.vel)
    if isa(wind.vel, Matrix) && size(wind.vel, 2) > 1
        wind_speed = wind.vel[1, 2]  # First time point, wind speed column
    elseif isa(wind.vel, Real)
        wind_speed = wind.vel
    else
        wind_speed = "Variable (see wind.vel)"
    end
    println("Free-flow wind speed: $wind_speed m/s")
else
    println("Free-flow wind speed: Not available")
end

if hasfield(typeof(wind), :ti) && !isnothing(wind.ti)
    if isa(wind.ti, Matrix) && size(wind.ti, 2) > 1
        turbulence_intensity = wind.ti[1, 2]  # First time point, TI column
    elseif isa(wind.ti, Real)
        turbulence_intensity = wind.ti
    else
        turbulence_intensity = "Variable (see wind.ti)"
    end
    println("Turbulence intensity: $(round(turbulence_intensity * 100, digits=1))%")
else
    println("Turbulence intensity: Not available")
end
