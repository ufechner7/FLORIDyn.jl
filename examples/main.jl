# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using Timers
tic()
using FLORIDyn, TerminalPager, ControlPlots

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
PARALLEL = false

# PLT options:
# PLT=1: Velocity reduction plot (if not using online visualization)
# PLT=2: Added turbulence plot  
# PLT=3: Wind speed plot
# PLT=4: Measurements plot (separated subplots)
# PLT=5: Measurements plot (combined)
# PLT=6: Velocity reduction plot with online visualization
# PLT=7: Create videos from saved frames
if !  @isdefined PLT; PLT=1; end

function get_parameters(vis)
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)

    # create settings struct
    set = Settings(wind, sim, con)
    set.parallel = PARALLEL

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

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions until no more change happens (wrong comment in original code)
wf = initSimulation(wf, sim)

toc()
# 0.115 s on Desktop, 0.39 s with MATLAB
# 0.115 seconds (891.24 k allocations: 368.147 MiB, 10.18% gc time)
# 0.110 seconds (883.14 k allocations: 272.819 MiB, 9.62% gc time) iterateOPs! allocation free
# 0.081 seconds (723.31 k allocations: 168.226 MiB, 8.10% gc time) findTurbineGroups allocation free

# @info "Type 'md |> pager' to see the results of the simulation."
# @info "Type 'q' to exit the pager."

if PLT == 1
    vis.online = false
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris)
    plotFlowField(plt, wf, X, Y, Z, vis; msr=1)
elseif PLT == 2
    vis.online = false
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris)
    plotFlowField(plt, wf, X, Y, Z, vis; msr=2)
elseif PLT == 3
    vis.online = false
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    @time Z, X, Y = calcFlowField(set, wf, wind, floris)
    plotFlowField(plt, wf, X, Y, Z, vis; msr=3)
elseif PLT == 4
    vis.online = false
    wf, md, set, floris, wind = get_parameters(vis)
    plotMeasurements(plt, wf, md, vis; separated=true)
elseif PLT == 5
    vis.online = false
    wf, md, set, floris, wind = get_parameters(vis)
    plotMeasurements(plt, wf, md, vis; separated=false)
elseif PLT == 6
    vis.online = true
    # Clean up any existing PNG files in video folder before starting
    if isdir("video")
        println("Cleaning up existing PNG files in video folder...")
        video_files = readdir("video")
        png_files = filter(f -> endswith(f, ".png"), video_files)
        for file in png_files
            try
                rm(joinpath("video", file))
            catch e
                @warn "Failed to delete $file: $e"
            end
        end
        if !isempty(png_files)
            println("Deleted $(length(png_files)) PNG files")
        end
    end
    @time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
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
