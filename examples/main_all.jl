# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Run all visualisations that are defined in vis_default.yaml
using Timers
using FLORIDyn, TerminalPager, DistributedNext, ControlPlots

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

vis = Vis(vis_file)

# Automatic parallel/threading setup
tic()
include("remote_plotting.jl")
toc()

FLORIDYN_EXECUTED = false

function get_parameters(vis)
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)

    # create settings struct with automatic parallel/threading detection
    set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)  
    wf = initSimulation(wf, sim)
    wf, md, mi = run_floridyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)
    return wf, md, set, floris, wind 
end

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

# Process flow fields
if length(vis.flow_fields) == 0
    @info "Skipping flow field visualisation."
end
for flow_field in vis.flow_fields
    global wf, md, mi, FLORIDYN_EXECUTED
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
                createAllVideos(fps=4, delete_frames=false, output_dir=vis.output_path)
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
    @info "Plotting measurements: $(measurement.name)"
    @time plot_measurements(wf, md, vis; separated=measurement.separated, plt)
end   


# if PLT == 1
#     vis.online = false
#     @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
#     @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
#     @time plot_flow_field(wf, X, Y, Z, vis; msr=1, plt)
# elseif PLT == 2
#     vis.online = false
#     @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
#     @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
#     @time plot_flow_field(wf, X, Y, Z, vis; msr=2, plt)
# elseif PLT == 3
#     vis.online = false
#     @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
#     @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
#     @time plot_flow_field(wf, X, Y, Z, vis; msr=3, plt)
# elseif PLT == 4
#     vis.online = false
#     wf, md, set, floris, wind = get_parameters(vis)
#     @time plot_measurements(wf, md, vis; separated=true, plt)
# elseif PLT == 5
#     vis.online = false
#     wf, md, set, floris, wind = get_parameters(vis)
#     @time plot_measurements(wf, md, vis; separated=false, plt)
# elseif PLT == 6
#     vis.online = true
#     # Clean up any existing PNG files in video folder before starting
#     cleanup_video_folder()
#     @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
# elseif PLT == 7
#     # Create videos from saved plot frames
#     println("Creating videos from saved plot frames...")
#     if isdir("video")
#         video_paths = createAllVideos(fps=4, delete_frames=false, output_dir=vis.output_path)
#         if !isempty(video_paths)
#             println("Videos created successfully!")
#             for path in video_paths
#                 println("  - $path")
#             end
#         else
#             println("No videos created. Make sure you have run the simulation with vis.save=true first.")
#         end
#     else
#         println("No 'video' directory found. Run simulation with vis.save=true to generate frames first.")
#     end
# end
nothing
