# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using Timers
tic()
using FLORIDyn, TerminalPager, Distributed

# PLT options:
# PLT=1: Velocity reduction plot
# PLT=2: Added turbulence plot  
# PLT=3: Wind speed plot
# PLT=4: Measurements plot (separated subplots)
# PLT=5: Measurements plot (combined)
# PLT=6: Velocity reduction plot with online visualization
# PLT=7: Create videos from saved frames
if !  @isdefined PLT; PLT=1; end
if !  @isdefined LAST_PLT; LAST_PLT=PLT; end

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
PARALLEL  = false
THREADING = false

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
    @everywhere function plot_flow_field(wf, X, Y, Z, vis; msr=3)
        # Create a fresh plt instance just for this task
        local_plt = ControlPlots.plt
        return plotFlowField(local_plt, wf, X, Y, Z, vis; msr=msr)
    end
    @everywhere function plot_measurements(wf, md, vis; separated)
        # Create a fresh plt instance just for this task
        local_plt = ControlPlots.plt
        return plotMeasurements(local_plt, wf, md, vis; separated=separated)
    end
    @everywhere function close_all()
        # Create a fresh plt instance just for this task
        local_plt = ControlPlots.plt
        return local_plt.close("all")
    end
    toc()
else
    @eval using ControlPlots
end

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

if LAST_PLT == PLT || LAST_PLT == 1 && PLT == 6 || LAST_PLT == 6 && PLT == 1
    if set.parallel
        @spawnat 2 close_all()
    else
        plt.close("all")
    end
end

if PLT == 1
    vis.online = false
    GC.gc()
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
    if set.parallel
        @time @spawnat 2 plot_flow_field(wf, X, Y, Z, vis; msr=1)
    else
        @time plotFlowField(plt, wf, X, Y, Z, vis; msr=1)
    end
elseif PLT == 2
    vis.online = false
    GC.gc()
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
    if set.parallel
        @time @spawnat 2 plot_flow_field(wf, X, Y, Z, vis; msr=2)
    else
        @time plotFlowField(plt, wf, X, Y, Z, vis; msr=2)
    end
elseif PLT == 3
    vis.online = false
    GC.gc()
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
        if set.parallel
        @time @spawnat 2 plot_flow_field(wf, X, Y, Z, vis; msr=3)
    else
        @time plotFlowField(plt, wf, X, Y, Z, vis; msr=3)
    end
elseif PLT == 4
    vis.online = false
    GC.gc()
    wf, md, set, floris, wind = get_parameters(vis)
    if set.parallel
        @time @spawnat 2 plot_measurements(wf, md, vis; separated=true)
    else
        plotMeasurements(plt, wf, md, vis; separated=true)
    end
elseif PLT == 5
    vis.online = false
    GC.gc()
    wf, md, set, floris, wind = get_parameters(vis)
    if set.parallel
        @time @spawnat 2 plot_measurements(wf, md, vis; separated=false)
    else
        plotMeasurements(plt, wf, md, vis; separated=false)
    end
elseif PLT == 6
    GC.gc()
    vis.online = true
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
LAST_PLT = PLT
nothing
