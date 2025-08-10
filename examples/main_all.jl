# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Run all visualisations that are defined in vis_default.yaml
using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext, ControlPlots, JLD2

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

vis = Vis(vis_file)

# Automatic parallel/threading setup
include("remote_plotting.jl")

FLORIDYN_EXECUTED = false

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions
wf = initSimulation(wf, sim)
@info "Initial conditions done, starting simulation..."
toc()

close_all(plt)
if vis.unique_output_folder
    vis.unique_folder = unique_name()
else
    vis.unique_folder = ""
end

# Process flow fields
if length(vis.flow_fields) == 0
    @info "Skipping flow field visualisation."
end
for flow_field in vis.flow_fields
    global wf, md, mi, FLORIDYN_EXECUTED
    
    # Skip this flow field if skip flag is set
    if flow_field.skip
        @info "Skipping flow field: $(flow_field.name)"
        continue
    end
    
    Z, X, Y = calcFlowField(set, wf, wind, floris)
    if flow_field.name == "flow_field_vel_reduction"
        msr = 1
    elseif flow_field.name == "flow_field_added_turbulence"
        msr = 2
    elseif flow_field.name == "flow_field_eff_wind_speed"
        msr = 3
    else
        msr = 1  # default to velocity reduction
    end
    if flow_field.online
        cleanup_video_folder()
        vis.online = true
        @info "Starting simulation with online visualisation for flow field $(flow_field.name) ..."
        wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
        FLORIDYN_EXECUTED = true

        if flow_field.create_video
            @info "Creating video for flow field $(flow_field.name) ..."
            video_paths = redirect_stdout(devnull) do
                createAllVideos(fps=4, delete_frames=false, video_dir=vis.video_path, output_dir=vis.output_path)
            end
            if !isempty(video_paths)
                @info "Video created successfully: $(video_paths[1])"
            else
                @warn "No video created."
            end
        end
    else
        vis.online = false
        @info "Plotting flow field: $(flow_field.name)"
        @time plot_flow_field(wf, X, Y, Z, vis; msr, plt)
    end
end
if length(vis.measurements) > 0 && ! FLORIDYN_EXECUTED
    global wf, md, mi
    vis.online = false
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)    
end
for measurement in vis.measurements
    # Skip this measurement if skip flag is set
    if measurement.skip
        @info "Skipping measurement: $(measurement.name), separated: $(measurement.separated)"
        continue
    end
    
    @info "Plotting measurement: $(measurement.name), separated: $(measurement.separated)"
    if measurement.name == "msr_vel_reduction"
        msr = 1
    elseif measurement.name == "msr_added_turbulence"
        msr = 2
    elseif measurement.name == "msr_eff_wind_speed"
        msr = 3
    else
        msr = 1  # default to velocity reduction
    end
    plot_measurements(wf, md, vis; separated=measurement.separated, msr, plt)
end

if vis.save_results
    # Save simulation results as .jld2 files
    @info "Saving simulation results..."
    results_filename = joinpath(vis.output_path, "results.jld2")
    
    try
        # Save the main simulation variables
        jldsave(results_filename; wf, md, mi)
        if vis.print_filenames
            @info "Simulation results saved to: $(results_filename)"
        end
    catch e
        @error "Failed to save simulation results: $e"
    end
end

toc()

nothing
