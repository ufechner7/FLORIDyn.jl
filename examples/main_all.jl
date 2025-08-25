# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Run all visualisations that are defined in vis_default.yaml
using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext, ControlPlots, JLD2
using Logging, LoggingExtras, Dates

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

vis = Vis(vis_file)

if (@isdefined plt) && !isnothing(plt)
    if !vis.show_plots
        plt.ioff()
        plt.pygui(false)
    else
        plt.ion()
    end
end

# Automatic parallel/threading setup
include("remote_plotting.jl")

FLORIDYN_EXECUTED = false

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

if vis.unique_output_folder
    vis.unique_folder = unique_name()
else
    vis.unique_folder = ""
end
log_path = joinpath(vis.output_path, "simulation.log")
base_file_logger = LoggingExtras.FormatLogger(log_path) do io, meta
    println(io, meta.level, ": ", meta.message)
end
# Keep only Info and above (drop Debug) unless explicitly enabled via vis.log_debug
min_level = vis.log_debug ? Logging.Debug : Logging.Info
file_logger = LoggingExtras.MinLevelLogger(base_file_logger, min_level)
console_logger = current_logger()
tee_logger = LoggingExtras.TeeLogger(console_logger, file_logger)
with_logger(tee_logger) do
    # Run initial conditions
    wf = initSimulation(wf, sim)
    @info "Initialization done after $(round(toc(false), digits=2)) s, starting simulation..."
    @warn "This may take a while, please be patient."

    close_all(plt)

    # Process flow fields
    if length(vis.flow_fields) == 0
        @info "Skipping flow field visualisation."
    end
    if ! vis.show_plots
        @info "Do not show the plots, only create png files."
    end
    for flow_field in vis.flow_fields
        global wf, md, mi, FLORIDYN_EXECUTED
        if vis.skip_flow_fields
            @info "Skipping flow field visualisation."
            break
        end
        if flow_field.skip
            @info "Skipping flow field: $(flow_field.name)"
            continue
        end
        Z, X, Y = calcFlowField(set, wf, wind, floris; vis)
        msr = toMSR(flow_field.name)
        if flow_field.online
            if ! vis.unique_output_folder
                @info "Cleaning up video folder: $(vis.video_path)"
                cleanup_video_folder()
            end
            vis.online = true
            @info "Starting simulation with online visualisation for $(flow_field.name) ..."
            wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris; msr)
            FLORIDYN_EXECUTED = true
            if flow_field.create_video
                @info "Creating video for $(flow_field.name) ..."
                video_paths = redirect_stdout(devnull) do
                    createAllVideos(fps=4, delete_frames=false, video_dir=vis.video_path, output_dir=vis.output_path)
                end
                if !isempty(video_paths)
                    @info "Video created successfully: $(video_paths[1])"
                    vis.no_videos += 1
                else
                    @warn "No video created."
                end
            end
        else
            vis.online = false
            @info "Plotting flow field: $(flow_field.name)"
            plot_flow_field(wf, X, Y, Z, vis; msr, plt)
            vis.no_plots += 1
        end
    end
    if length(vis.measurements) > 0 && ! FLORIDYN_EXECUTED
        global wf, md, mi
        vis.online = false
        wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    end
    for measurement in vis.measurements
        if measurement.skip
            @info "Skipping measurement: $(measurement.name), separated: $(measurement.separated)"
            continue
        end
        if vis.skip_measurements
            @info "Skipping measurement visualisation."
            break
        end
        @info "Plotting measurement: $(measurement.name), separated: $(measurement.separated)"
        msr = measurement.name == "msr_vel_reduction" ? VelReduction : measurement.name == "msr_added_turbulence" ? 
                                AddedTurbulence : measurement.name == "msr_eff_wind_speed" ? EffWind : VelReduction
        plot_measurements(wf, md, vis; separated=measurement.separated, msr, plt)
        vis.no_plots += 1
    end
    if vis.save_results
        @info "Saving simulation results..."
        results_filename = joinpath(vis.output_path, "results.jld2")
        try
            jldsave(results_filename; wf, md, mi)
            if vis.print_filenames
                @info "Simulation results saved to: $(results_filename)"
            end
        catch e
            @error "Failed to save simulation results: $e"
        end
    end
    shorten(p) = replace(p, homedir() => "~")
    target_dir = vis.output_path
    @info "Copying yaml files: $(shorten(settings_file)) -> $(shorten(target_dir))"
    @info "Copying yaml files: $(shorten(vis_file)) -> $(shorten(target_dir))"
    cp(settings_file, joinpath(target_dir, basename(settings_file)))
    cp(vis_file,      joinpath(target_dir, basename(vis_file)))
    try
        settings_basename = basename(settings_file)
        settings_foldername = first(splitext(settings_basename))
        src_dir = joinpath(dirname(settings_file), settings_foldername)
        dest_dir = joinpath(vis.output_path, settings_foldername)
        if isdir(src_dir)
            @info "Copying CSV folder: $(shorten(src_dir)) -> $(shorten(dest_dir))"
            cp(src_dir, dest_dir; force=true)
        else
            @info "No matching CSV folder to copy (expected directory not found): $(shorten(src_dir))"
        end
    catch e
        @warn "Failed to copy settings data folder" exception=(e, catch_backtrace())
    end
    @info "Simulation completed. Plots created: $(vis.no_plots), videos created: $(vis.no_videos) ."
    @info "Total execution time: $(round(toc(false), digits=2)) s"
end
nothing
