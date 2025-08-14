# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause
 
using FLORIDyn, Test, Statistics

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

 @testset "setUpTmpWFAndRun_basic" begin
    settings_file = "data/2021_9T_Data.yaml"
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    # create settings struct
    set = Settings(wind, sim, con)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    wf.dep = findTurbineGroups(wf, floridyn)
    wf.intOPs = interpolateOPs(wf)
    wf_old = deepcopy(wf)
    M, wf = setUpTmpWFAndRun(set, wf, floris, wind)
    @test ! structs_equal(wf_old, wf; prn=false)
end
nothing