# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Minimal example of how to run a simulation using FLORIDyn.jl 
# for benchmarking the 54 turbine layout.
using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext, Statistics 
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"

USE_STEP = true
USE_FEED_FORWARD = true
USE_MPC = false  # If false, use simple step control
ONLINE = false
PLOT_STEP_RESPONSE = false
PLOT_STORAGE_VS_WINDDIR = true
WIND_DIR = 270.0  # Wind direction for step response simulation
if PLOT_STEP_RESPONSE
    WIND_DIRS = 200:10:340  # Wind directions for step response simulation
else
    WIND_DIRS = 200:1:340  # Wind directions for storage vs wind direction simulation
end

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

function calc_demand_and_power(settings_file; wind_dir=WIND_DIR)
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    sim.end_time += 420
    con.yaw="Constant"
    wind.input_dir="Constant"
    wind.dir_fixed = wind_dir
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

    vis.online = ONLINE
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

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
    return rel_power, demand_values, times, wind, sim
end

function analyse_results(rel_power, demand_values; dt=sim.time_step)
    t1 = 600
    t2 = 700
    t3 = 1200
    t4 = 1600
    mean_peak = mean(rel_power[1+t1÷dt:1+t2÷dt])
    mean_final = mean(rel_power[1+t3÷dt:1+t4÷dt])
    extra_power = mean_peak - mean_final
    println("\n--- Analysis of Power Output ---")
    println("Relative peak power:  $(round(mean_peak * 100, digits=2))%")
    println("Relative final power: $(round(mean_final * 100, digits=2))%")
    println("Extra power:          $(round(extra_power * 100, digits=2))%")
    storage_time = sum((rel_power[1+t1÷dt:1+t4÷dt] .- mean_final) .* dt)
    println("Storage time at full power: $(round(storage_time, digits=2))s")
end

function step_response(wind_dirs=WIND_DIRS)
    # First, calculate for the single WIND_DIR for detailed analysis
    rel_power, demand_values, times, wind, sim = calc_demand_and_power(settings_file; wind_dir=WIND_DIR)

    plot_rmt(times, [rel_power .* 100, demand_values .* 100]; xlabel="Time [s]", 
             xlims=(400+sim.time_step, 1200+sim.time_step), ylabel="Rel. Power Output [%]", 
             labels=["rel_power", "rel_demand"], fig="Step Response - Wind Dir $(WIND_DIR)°", pltctrl=pltctrl)

    # Save the first plot
    if pltctrl !== nothing
        filename1 = "docs/src/step_response_wind_dir_$(WIND_DIR).png"
        mkpath(dirname(filename1))
        try
            pltctrl.plt.savefig(filename1, dpi=150, bbox_inches="tight", pad_inches=0.1, facecolor="white")
            println("Saved plot: $filename1")
        catch e
            @warn "Failed to save plot: $e"
        end
    else
        @warn "Saving the plot only works in single-threaded mode, launch Julia with jl2!"
    end

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

    analyse_results(rel_power, demand_values; dt=sim.time_step)

    # Now calculate and plot relative power for all wind directions
    println("\n--- Calculating relative power for all wind directions ---")
    
    # Storage for results  
    all_rel_powers = Vector{Float64}[]  # Properly typed as Vector of Vector{Float64}
    labels = String[]
    
    # Add demand curve for reference
    push!(all_rel_powers, demand_values .* 100)
    push!(labels, "demand")
    
    # Calculate for each wind direction
    for wd in wind_dirs
        println("Processing wind direction: $(wd)°")
        rel_power_wd, _, _, _, _ = calc_demand_and_power(settings_file; wind_dir=wd)
        push!(all_rel_powers, rel_power_wd .* 100)
        push!(labels, "$(wd)°")
    end
    
    # Plot all relative power curves
    plot_rmt(times, all_rel_powers; 
             xlabel="Time [s]", 
             xlims=(400+sim.time_step, 1200+sim.time_step),
             ylabel="Rel. Power Output [%]", 
             labels=labels, 
             fig="Step Response - All Wind Directions",
             pltctrl=pltctrl)

    # Save the second plot
    if pltctrl !== nothing
        filename2 = "docs/src/step_response_all_wind_directions.png"
        mkpath(dirname(filename2))
        try
            pltctrl.plt.savefig(filename2, dpi=150, bbox_inches="tight", pad_inches=0.1, facecolor="white")
            println("Saved plot: $filename2")
        catch e
            @warn "Failed to save plot: $e"
        end
    end
             
    println("Completed step response analysis for all wind directions: $(collect(wind_dirs))°")
    
    # Execute the documentation script
    println("\n--- Generating documentation ---")
    try
        # Change to the project root directory and execute the documentation script
        script_path = joinpath(pwd(), "bin", "document_examples")
        if isfile(script_path)
            run(`bash $script_path`)
            println("Documentation generated successfully")
        else
            @warn "Documentation script not found at: $script_path"
        end
    catch e
        @warn "Failed to execute documentation script: $e"
    end
end

function storage_vs_winddir(settings_file; wind_dirs= WIND_DIRS)
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    dt = sim.time_step  # seconds
    extra_powers = Float64[]
    storage_times = Float64[]
    for wd in wind_dirs
            println("\n--- Wind Direction: $(wd) ° ---")
            rel_power, demand_values, times, wind, sim = calc_demand_and_power(settings_file; wind_dir=wd)
        t1 = 600
        t2 = 700
        t3 = 1200
        t4 = 1600
        mean_peak = mean(rel_power[1+t1÷sim.time_step:1+t2÷sim.time_step])
        @assert mean_peak > 0.98 "Mean peak power < 98%"
        mean_final = mean(rel_power[1+t3÷sim.time_step:1+t4÷sim.time_step])
        extra_power = mean_peak - mean_final
        storage_time = sum((rel_power[1+t1÷dt:1+t4÷dt] .- mean_final) .* dt)
        push!(extra_powers, extra_power)
        push!(storage_times, storage_time)
        println("Extra power: $(round(extra_power * 100, digits=2))%")
        println("Storage time at full power: $(round(storage_time, digits=2))s")
    end
    plot_rmt(wind_dirs, extra_powers .* 100; xlabel="Wind Direction [°]", ylabel="Extra Power [%]", fig="Extra Power", pltctrl=pltctrl)
    
    # Save the extra power plot
    if pltctrl !== nothing
        filename_extra = "docs/src/extra_power_vs_wind_dir.png"
        mkpath(dirname(filename_extra))
        try
            pltctrl.plt.savefig(filename_extra, dpi=150, bbox_inches="tight", pad_inches=0.1, facecolor="white")
            println("Saved plot: $filename_extra")
        catch e
            @warn "Failed to save extra power plot: $e"
        end
    end
    
    plot_rmt(wind_dirs, storage_times; xlabel="Wind Direction [°]", ylabel="Storage Time at Full Power [s]", fig="Storage Time", pltctrl=pltctrl)
    
    # Save the storage time plot
    if pltctrl !== nothing
        filename_storage = "docs/src/storage_time_vs_wind_dir.png"
        mkpath(dirname(filename_storage))
        try
            pltctrl.plt.savefig(filename_storage, dpi=150, bbox_inches="tight", pad_inches=0.1, facecolor="white")
            println("Saved plot: $filename_storage")
        catch e
            @warn "Failed to save storage time plot: $e"
        end
    else
        @warn "Saving the plot only works in single-threaded mode, launch Julia with jl2!"
    end
end

# Run the step response simulation
if PLOT_STEP_RESPONSE
    step_response()
end
# Run the storage vs wind direction simulation
if PLOT_STORAGE_VS_WINDDIR
    storage_vs_winddir(settings_file)
end