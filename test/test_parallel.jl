# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

if !isdefined(Main, :Test)
    using Test
end 

if !isdefined(Main, :DistributedNext)
    using DistributedNext
end 

if ! isinteractive()
if !isdefined(Main, :FLORIDyn)
    using FLORIDyn
end

if !isdefined(Main, :ControlPlots)
    using ControlPlots
    @info "using ControlPlots"
end

function get_parameters(vis, settings_file)
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)

    # create settings struct with automatic parallel/threading detection
    set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)  
    wf = initSimulation(wf, sim)
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    return wf, md, set, floris, wind 
end
@testset "multithreading" begin
    @test FLORIDyn.nthreads() > 1
    settings_file = "data/2021_9T_Data.yaml"
    vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
    vis.online = false
    vis.unit_test = true
    for i in 1:8
        local wf, md, set, floris, wind, X, Y, Z
        wf, md, set, floris, wind = get_parameters(vis, settings_file)
        @time Z, X, Y = calcFlowField(set, wf, wind, floris; plt)
        msr = mod(i - 1, 3) + 1  # Convert to 1-based indexing (1, 2, 3, 1, 2, 3)
        plot_flow_field(wf, X, Y, Z, vis; msr, plt=ControlPlots.plt)
        @test true
        GC.gc()  # Force garbage collection between iterations
    end
end
@testset "parallel" begin
    # this is now done in runtests.jl
    # include("../src/visualisation/remote_plotting.jl") 
    init_plotting()
    @test workers()[1] >= 2
end
else
    # Running tests via Pkg.test (safest approach)
    @eval Main using Pkg
    Main.Pkg.test(test_args=["test_parallel.jl"])
end
nothing