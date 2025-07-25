# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test
@testset verbose=true "floridyncl" begin
    include("test_prepare_simulation.jl")
    @testset "angSOWFA2world" begin
        # Test 1: deg_SOWFA = 0 should give rad_World = deg2rad(270)
        @test isapprox(angSOWFA2world(0), deg2rad(270))

        # Test 2: deg_SOWFA = 90 should give rad_World = deg2rad(180)
        @test isapprox(angSOWFA2world(90), deg2rad(180))

        # Test 3: deg_SOWFA = 180 should give rad_World = deg2rad(90)
        @test isapprox(angSOWFA2world(180), deg2rad(90))

        # Test 4: deg_SOWFA = 270 should give rad_World = deg2rad(0)
        @test isapprox(angSOWFA2world(270), deg2rad(0))

        # Test 5: deg_SOWFA = 360 should give rad_World = deg2rad(-90)
        @test isapprox(angSOWFA2world(360), deg2rad(-90))
    end
    function structs_equal(a::T, b::T; prn=true) where T
        result = true
        fields = fieldnames(T)
        for f in fields
            val_a = getfield(a, f)
            val_b = getfield(b, f)
            if val_a != val_b
                prn && println("Field $(f): a = $(val_a), b = $(val_b)")
                result = false
            end
        end
        return result
    end
    @testset "initSimulation" begin
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        # % Load linked data
        turbProp        = turbineArrayProperties(settings_file)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbProp, sim)
        wf_old = deepcopy(wf)
        wf = initSimulation(wf, sim)
        @test structs_equal(wf_old, wf)
        path = joinpath(sim.path_to_data, "T_init.jld2")
        rm(path; force=true)
        sim.init="init"
        sim.save_init_state = true
        wf = initSimulation(wf, sim)
        @test isfile(path)
        sim.init="load"
        sim.save_init_state = false
        wf = nothing
        wf = initSimulation(wf, sim)
        @test structs_equal(wf_old, wf)
        rm(path; force=true)
    end
    @testset "perturbationOfTheWF!" begin
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        # % Load linked data
        turbProp        = turbineArrayProperties(settings_file)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbProp, sim)
        wf_old = deepcopy(wf)
        perturbationOfTheWF!(wf, wind)
        @test structs_equal(wf_old, wf)
        wind.pertubation.vel = true
        perturbationOfTheWF!(wf, wind)
        @test ! structs_equal(wf_old, wf; prn=false)
        wind.pertubation.vel = false
        wind.pertubation.dir = true
        perturbationOfTheWF!(wf, wind)
        @test ! structs_equal(wf_old, wf; prn=false)
        wind.pertubation.dir = false
        wind.pertubation.ti = true
        perturbationOfTheWF!(wf, wind)
        @test ! structs_equal(wf_old, wf; prn=false)
    end
    @testset "findTurbineGroups" begin
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        # % Load linked data
        turbProp        = turbineArrayProperties(settings_file)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbProp, sim)
        vv_dep = findTurbineGroups(wf, floridyn)
        vv_dep_expected =   [
                                Int[],
                                [1],
                                [1, 2, 4],
                                Int[],
                                [4],
                                [4, 5, 7],
                                Int[],
                                [7],
                                [7, 8]
                            ]
        @test vv_dep == vv_dep_expected
    end
    @testset "iterateOPs!" begin
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        # % Load linked data
        turbine_prop        = turbineArrayProperties(settings_file)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_prop, sim)
        wf_old = deepcopy(wf)
        iterateOPs!(set.iterate_mode, wf, sim, floris, floridyn)
        @test ! structs_equal(wf_old, wf; prn=false)
    end
    @testset "interpolateOPs" begin
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        # % Load linked data
        turbine_prop        = turbineArrayProperties(settings_file)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_prop, sim)
        wf.dep = findTurbineGroups(wf, floridyn)
        wf.intOPs = interpolateOPs(wf)
        @test length(wf.intOPs) == wf.nT
    end
    @testset "setUpTmpWFAndRun" begin
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        # % Load linked data
        turbine_prop        = turbineArrayProperties(settings_file)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_prop, sim)
        wf.dep = findTurbineGroups(wf, floridyn)
        wf.intOPs = interpolateOPs(wf)
        wf_old = deepcopy(wf)
        M, wf = setUpTmpWFAndRun(set, wf, floris, wind)
        @test ! structs_equal(wf_old, wf; prn=false)
    end
    @testset "runFLORIDyn" begin
        global md
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        # % Load linked data
        turbine_prop        = turbineArrayProperties(settings_file)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_prop, sim)
        wf, md, mi = runFLORIDyn(set, wf, wind, sim, con, floridyn, floris)
        @test size(md) == (2709, 6) # from Matlab
        @test minimum(md.ForeignReduction) ≈ 72.57019949691814 # Matlab: 73.8438
        @test mean(md.ForeignReduction)    ≈ 98.54434468415639 # Matlab: 98.2902
    end
end # testset floridyncl
nothing