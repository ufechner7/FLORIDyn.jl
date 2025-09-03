# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    calc_rel_power(settings_file; dt=350, wind_dir=180.0)

Calculate relative power output for a wind farm simulation over time.

This function runs a FLORIDyn simulation and computes the relative power output based on 
the cube of relative wind speeds at each turbine. The relative power represents the 
power output compared to undisturbed wind conditions.

# Arguments
- `settings_file::String`: Path to the YAML configuration file containing wind farm settings
- `dt::Real=350`: Additional simulation time in seconds beyond the base configuration
- `wind_dir::Union{Real,Nothing}=180.0`: Wind direction in degrees. If `nothing`, uses 
  variable wind direction from the settings file. If a number, uses constant wind direction.

# Returns
A tuple containing:
- `times::Vector{Float64}`: Time points of the simulation [s]
- `rel_power::Vector{Float64}`: Relative power output at each time point [dimensionless]
- `set::Settings`: Simulation settings used
- `wf::WindFarm`: Wind farm object with final turbine states
- `wind::Wind`: Wind field object
- `floris::Floris`: FLORIS model object

# Physics
The relative power is calculated as:
```
rel_power[t] = (1/nT) * Σᵢ (rel_speed[i,t])³
```
where `rel_speed[i,t]` is the relative wind speed at turbine `i` and time `t`, and 
`nT` is the number of turbines. This follows the cubic relationship between wind 
speed and power: P ∝ v³.

# Threading
The function automatically detects and uses multithreading if available 
(`Threads.nthreads() > 1`). When using fixed wind direction, the simulation 
uses constant yaw and direction modes for improved performance.

# Example
```julia
# Variable wind direction simulation
times, rel_power, set, wf, wind, floris = calc_rel_power("data/2021_9T_Data.yaml"; 
                                                         dt=200, wind_dir=nothing)

# Fixed wind direction simulation  
times, rel_power, set, wf, wind, floris = calc_rel_power("data/2021_9T_Data.yaml"; 
                                                         dt=100, wind_dir=270.0)

# Analyze results
mean_power = mean(rel_power)
power_std = std(rel_power)
```

# See also
[`FLORIDyn.prepare_large_plot_inputs`](@ref), [`run_floridyn`](@ref), [`Settings`](@ref)
"""
function calc_rel_power(settings_file; dt=350, wind_dir=180.0, ti=0.062)
    fixed_wind_dir = ! isnothing(wind_dir)
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    sim.end_time += dt
    if fixed_wind_dir
        con.yaw = "Constant"
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
        con.yaw_data[1,1] = wind_dir
        wind.dir[1,1]     = wind_dir
    end
    if !isnothing(ti)
        wind.ti = ti
        # for iT in 1:wf.nT
        #     wf.ti[iT] = ti
        # end
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