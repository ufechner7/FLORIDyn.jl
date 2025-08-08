# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using Timers
using FLORIDyn, TerminalPager, DistributedNext, ControlPlots

# PLT options:
# PLT=1: Velocity reduction plot
# PLT=2: Added turbulence plot  
# PLT=3: Wind speed plot
# PLT=4: Measurements plot (separated subplots)
# PLT=5: Measurements plot (combined)
# PLT=6: Velocity reduction plot with online visualization
# PLT=7: Create videos from saved frames
if !  @isdefined PLT; PLT=1; end
if PLT == 6; NEW_PLT = 1; else NEW_PLT = PLT; end
if ! @isdefined LAST_PLT; LAST_PLT=Set(NEW_PLT); end

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"

vis = Vis(vis_file)

# Automatic parallel/threading setup
tic()
include("remote_plotting.jl")
toc()

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

# Run initial conditions until no more change happens (wrong comment in original code)
wf = initSimulation(wf, sim)
@info "Initial conditions done, starting simulation..."
toc()

if NEW_PLT in LAST_PLT
    # If the new plot was displayed before, close all plots
    close_all(plt)
end

if PLT == 1
    vis.online = false
    @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
    @time plot_flow_field(wf, X, Y, Z, vis; msr=1, plt)
elseif PLT == 2
    vis.online = false
    @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
    @time plot_flow_field(wf, X, Y, Z, vis; msr=2, plt)
elseif PLT == 3
    vis.online = false
    @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
    @time plot_flow_field(wf, X, Y, Z, vis; msr=3, plt)
elseif PLT == 4
    vis.online = false
    wf, md, set, floris, wind = get_parameters(vis)
    @time plot_measurements(wf, md, vis; separated=true, plt)
elseif PLT == 5
    vis.online = false
    wf, md, set, floris, wind = get_parameters(vis)
    @time plot_measurements(wf, md, vis; separated=false, plt)
elseif PLT == 6
    vis.online = true
    # Clean up any existing PNG files in video folder before starting
    cleanup_video_folder()
    @time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
elseif PLT == 7
    # Create videos from saved plot frames
    println("Creating videos from saved plot frames...")
    if isdir("video")
        video_paths = createAllVideos(fps=4, delete_frames=false, output_dir=vis.output_path)
        if !isempty(video_paths)
            println("Videos created successfully!")
            for path in video_paths
                println("  - $path")
            end
        else
            println("No videos created. Make sure you have run the simulation with vis.save=true first.")
        end
    else
        println("No 'video' directory found. Run simulation with vis.save=true to generate frames first.")
    end
end
push!(LAST_PLT, NEW_PLT)
nothing
