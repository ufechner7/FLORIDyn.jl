# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

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
    wf, md, mi = run_floridyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)

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
    return times, rel_power, set, wf, wind, floris
end