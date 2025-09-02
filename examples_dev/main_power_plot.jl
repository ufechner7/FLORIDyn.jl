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

function calc_rel_power(settings_file; dt=350, wind_dir=180.0)
    fixed_wind_dir = ! isnothing(wind_dir)
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    sim.end_time += dt
    if fixed_wind_dir
        con.yaw = "Constant"
        con.yaw_data = [wind_dir;;]
        wind.input_dir = "Constant"
    end

    # create settings struct with automatic parallel/threading detection
    set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
    if fixed_wind_dir
        set.dir_mode = Direction_Constant()
        set.control_mode = Yaw_Constant()
    end

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    if fixed_wind_dir
        wind.dir[1,1] = wind_dir
    end

    vis = Vis()
    vis.online = false
    @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

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
    return times, rel_power
end

times, rel_power = calc_rel_power(settings_file; dt=350, wind_dir=nothing)

p = plot_rmt(times, rel_power .* 100; xlabel="Time [s]", ylabel="Rel. Power Output [%]", pltctrl)
display(p)

println("\nMean Relative Power Output:  $(round((mean(rel_power) * 100), digits=2)) %")
println("Final Relative Power Output: $(round((rel_power[end] * 100), digits=2)) %")
