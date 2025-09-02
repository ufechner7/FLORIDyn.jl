# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example how to plot the average, relative wind park power
# Either a fixed or a variable wind direction can be used 
using FLORIDyn, TerminalPager, DistributedNext, Statistics, JLD2
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"
WIND_DIR_MIN        = 270-90
WIND_DIR_MAX        = 270+90
WIND_DIR_STEPS      = Int(180/2.5) # 2.5° steps
LOAD_RESULTS        = false
SAVE_RESULTS        = true

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

function calc_pwr(wind_dir)
    times, rel_power, set, wf, wind, floris = calc_rel_power(settings_file; dt=350, wind_dir=wind_dir)
    return mean(rel_power), rel_power[end]
end

if WIND_DIR_MIN == WIND_DIR_MAX
    global wind_dirs, mean_pwrs, final_pwrs
    local mean_pwrs, final_pwrs
    mean_pwr, final_pwr = calc_pwr(WIND_DIR_MIN)
    wind_dirs = [WIND_DIR_MIN]
    mean_pwrs = [mean_pwr]
    final_pwrs = [final_pwr]
else
    global wind_dirs, mean_pwrs, final_pwrs
    local mean_pwrs, final_pwrs
    wind_dirs = collect(LinRange(WIND_DIR_MIN, WIND_DIR_MAX, WIND_DIR_STEPS))
    mean_pwrs = Float64[]
    final_pwrs = Float64[]
    for wd in wind_dirs
        mean_pwr, final_pwr = calc_pwr(wd)
        push!(mean_pwrs, mean_pwr)
        push!(final_pwrs, final_pwr)
    end
end

if SAVE_RESULTS
    @info "Saving simulation results..."
    results_filename = joinpath(vis.output_path, "results_power_vs_wind_dir.jld2")
    try
        jldsave(results_filename; wind_dirs, mean_pwrs, final_pwrs)
        if vis.print_filenames
            @info "Simulation results saved to: $(results_filename)"
        end
    catch e
        @error "Failed to save simulation results: $e"
    end
end

plot_rmt(wind_dirs, final_pwrs; xlabel="Wind Direction (deg)", ylabel="Relative Power", 
         title="Relative Windfarm Power vs Wind Direction")


