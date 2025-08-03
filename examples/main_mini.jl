# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using Timers
tic()
using FLORIDyn, TerminalPager, DistributedNext, ControlPlots

# PLT options:
# PLT=1: Velocity reduction plot
if !  @isdefined PLT; PLT=1; end
if PLT == 6; NEW_PLT = 1; else NEW_PLT = PLT; end
if ! @isdefined LAST_PLT; LAST_PLT=Set(NEW_PLT); end

settings_file = "data/2021_9T_Data.yaml"
vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)

# Automatic parallel/threading setup
if Threads.nthreads() > 1
    tic()
    include("../src/visualisation/remote_plotting.jl") 
    init_plotting()  # This sets up workers and remote plotting capabilities
    # Create a completely isolated plt instance for this specific task
    @everywhere function plot_flow_field(wf, X, Y, Z, vis; msr=3)
        # Create a fresh plt instance just for this task
        local_plt = ControlPlots.plt
        return plotFlowField(local_plt, wf, X, Y, Z, vis; msr=msr)
    end
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

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions until no more change happens (wrong comment in original code)
wf = initSimulation(wf, sim)
toc()

vis.online = false
GC.gc()
@time wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
@time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
@time smart_plot_flow_field(wf, X, Y, Z, vis; msr=1, plt=plt)

push!(LAST_PLT, NEW_PLT)
nothing
