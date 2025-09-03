# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example how to plot the average, relative wind park power
# Either a fixed or a variable wind direction can be used 
using FLORIDyn, TerminalPager, DistributedNext, Statistics
import JLD2
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"
TI             = 0.062*1  # turbulence intensity for fixed wind direction simulations
TIs            = [0.0, 0.062, 0.062*2]
WIND_DIR_MIN        = 270-90
WIND_DIR_MAX        = 270+90
WIND_DIR_STEPS      = Int(180/2.5)+1
LOAD_RESULTS        = false   # set to false to always run new simulations
RUN_SIMULATION      = false   # set to true to run simulations
SAVE_RESULTS        = false
PLOT_TIs            = true  # set to false to avoid loading non-existent TI results

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

# Initialize global variables
wind_dirs = Float64[]
mean_pwrs = Float64[]
final_pwrs = Float64[]

results_filename = joinpath(vis.output_path, "results_power_vs_wind_dir_ti_$(TI).jld2")

if LOAD_RESULTS
    @info "Loading previously saved simulation results..."
    try
        global wind_dirs, mean_pwrs, final_pwrs, RUN_SIMULATION, TI
        RUN_SIMULATION = false
        wind_dirs = JLD2.load(results_filename, "wind_dirs")
        mean_pwrs = JLD2.load(results_filename, "mean_pwrs")
        final_pwrs = JLD2.load(results_filename, "final_pwrs")
        TI = JLD2.load(results_filename, "TI")
    catch e
        @warn "Failed to load simulation results: $e"
        @info "Proceeding to run new simulations..."
        # Initialize empty arrays if loading fails
        global wind_dirs, mean_pwrs, final_pwrs, RUN_SIMULATION
        RUN_SIMULATION = true
        wind_dirs = Float64[]
        mean_pwrs = Float64[]
        final_pwrs = Float64[]
    end
end
if RUN_SIMULATION
    function calc_pwr(wind_dir)
        times, rel_power, set, wf, wind, floris = calc_rel_power(settings_file; dt=350, wind_dir=wind_dir)
        return mean(rel_power), rel_power[end]
    end

    global wind_dirs, mean_pwrs, final_pwrs
    if WIND_DIR_MIN == WIND_DIR_MAX
        mean_pwr, final_pwr = calc_pwr(WIND_DIR_MIN)
        wind_dirs = [WIND_DIR_MIN]
        mean_pwrs = [mean_pwr]
        final_pwrs = [final_pwr]
    else
        wind_dirs = collect(LinRange(WIND_DIR_MIN, WIND_DIR_MAX, WIND_DIR_STEPS))
        mean_pwrs = Float64[]
        final_pwrs = Float64[]
        for wd in wind_dirs
            @info "wind_dir: $(wd)"
            mean_pwr_, final_pwr_ = calc_pwr(wd)
            push!(mean_pwrs, mean_pwr_)
            push!(final_pwrs, final_pwr_)
        end
    end
end

if SAVE_RESULTS && RUN_SIMULATION # save the results only if a simulation run happened
    @info "Saving simulation results..."
    try
        jldsave(results_filename; wind_dirs, mean_pwrs, final_pwrs, TI)
        if vis.print_filenames
            @info "Simulation results saved to: $(results_filename)"
        end
    catch e
        @error "Failed to save simulation results: $e"
    end
end

wind_dirs_ = Float64[]  # Initialize as global 
final_pwr_vec = Vector{Float64}[]
labels = Vector{String}()
if PLOT_TIs
    global wind_dirs_, final_pwr_vec, labels  # Declare as global to avoid scoping issues
    for ti in TIs
        local results_filename_
        @info "Plotting results for TI=$(100*ti) %..."
        results_filename_ = joinpath(vis.output_path, "results_power_vs_wind_dir_ti_$(ti).jld2")
        println("Loading from $results_filename_ ...")
        try
            global wind_dirs_, final_pwr_vec  # Ensure these are treated as global
            wind_dirs_  = JLD2.load(results_filename_, "wind_dirs")
            final_pwrs_ = JLD2.load(results_filename_, "final_pwrs")
            TI_         = JLD2.load(results_filename_, "TI")
            @assert TI_ == ti
            push!(final_pwr_vec, final_pwrs_)
            push!(labels, "TI: $(100*ti) %")
        catch e
            @warn "Failed to load simulation results for TI=$(100*ti) %: $e"
            continue
        end
    end
    if !isempty(wind_dirs_)
        plot_rmt(wind_dirs_, final_pwr_vec; xlabel="Wind Direction (deg)", ylabel="Relative Power", 
                 labels, fig="Mean Relative Windfarm Power vs Wind Direction", pltctrl=pltctrl)
    else
        @warn "No data to plot - all TI result files failed to load"
    end
else
    @info "Plotting results for TI=$(100*TI) %..."
    plot_rmt(wind_dirs, mean_pwrs; xlabel="Wind Direction (deg)", ylabel="Relative Power", xlims=(WIND_DIR_MIN, WIND_DIR_MAX),
             title="Mean Relative Windfarm Power vs Wind Direction", fig="TI: $(100*TI) %", pltctrl=pltctrl)
end

