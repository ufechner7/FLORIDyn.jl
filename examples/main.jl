# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using Timers
tic()
using FLORIDyn, TerminalPager, ControlPlots, Distributed

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
PARALLEL = true
THREADING = true

if PARALLEL
    tic()
    include("../src/visualisation/remote_plotting.jl") 
    init_plotting()  # This now returns the main process plt and creates plt on workers
    @everywhere function plot_flow_field(wf, X, Y, Z, vis, t_rel; msr=1)
        global plot_state
        if abs(t_rel) < 1e-6
            plot_state = nothing
        end
        local_plt = ControlPlots.plt
        plot_state = plotFlowField(plot_state, local_plt, wf, X, Y, Z, vis, t_rel; msr=msr)
        nothing
    end
    toc()
end

# PLT options:
# PLT=1: Velocity reduction plot (if not using online visualization)
# PLT=2: Added turbulence plot  
# PLT=3: Wind speed plot
# PLT=4: Measurements plot (separated subplots)
# PLT=5: Measurements plot (combined)
# PLT=6: Velocity reduction plot with online visualization
# PLT=7: Create videos from saved frames
if !  @isdefined PLT; PLT=1; end

function get_parameters(vis, parallel=PARALLEL, threading=THREADING)
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)

    # create settings struct
    set = Settings(wind, sim, con)
    set.parallel = parallel
    set.threading = threading

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)  
    wf = initSimulation(wf, sim)
    wf, md, mi = runFLORIDyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)
    return wf, md, set, floris, wind 
end

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct
set = Settings(wind, sim, con)
set.parallel = PARALLEL
set.threading = THREADING

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions until no more change happens (wrong comment in original code)
wf = initSimulation(wf, sim)
toc()

if PLT == 1
    vis.online = false
    set.parallel = false  # Disable parallel plotting for this case
    GC.gc()
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
    @time plotFlowField(plt, wf, X, Y, Z, vis; msr=1)
elseif PLT == 2
    vis.online = false
    set.parallel = false  # Disable parallel plotting for this case
    GC.gc()
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
    plotFlowField(plt, wf, X, Y, Z, vis; msr=2)
elseif PLT == 3
    vis.online = false
    set.parallel = false  # Disable parallel plotting for this case
    GC.gc()
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
    plotFlowField(plt, wf, X, Y, Z, vis; msr=3)
elseif PLT == 4
    vis.online = false
    set.parallel = false  # Disable parallel plotting for this case
    wf, md, set, floris, wind = get_parameters(vis)
    plotMeasurements(plt, wf, md, vis; separated=true)
elseif PLT == 5
    vis.online = false
    set.parallel = false  # Disable parallel plotting for this case
    wf, md, set, floris, wind = get_parameters(vis)
    plotMeasurements(plt, wf, md, vis; separated=false)
elseif PLT == 6
    vis.online = true
    set.parallel = true  # Enable parallel plotting for this case
    # Clean up any existing PNG files in video folder before starting
    cleanup_video_folder()
    if PARALLEL
        @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris, plot_flow_field)
    else
        @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    end
elseif PLT == 7
    # Create videos from saved plot frames
    println("Creating videos from saved plot frames...")
    if isdir("video")
        video_paths = createAllVideos(fps=4, delete_frames=false)
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
nothing
