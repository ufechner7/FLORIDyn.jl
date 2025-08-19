# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test

@testset verbose=true "toMSR" begin
    # Enum names
    @test FLORIDyn.toMSR("VelReduction") === FLORIDyn.VelReduction
    @test FLORIDyn.toMSR("AddedTurbulence") === FLORIDyn.AddedTurbulence
    @test FLORIDyn.toMSR("EffWind") === FLORIDyn.EffWind

    # Flow field names
    @test FLORIDyn.toMSR("flow_field_vel_reduction") === FLORIDyn.VelReduction
    @test FLORIDyn.toMSR("flow_field_added_turbulence") === FLORIDyn.AddedTurbulence
    @test FLORIDyn.toMSR("flow_field_eff_wind_speed") === FLORIDyn.EffWind

    # Measurement names
    @test FLORIDyn.toMSR("msr_vel_reduction") === FLORIDyn.VelReduction
    @test FLORIDyn.toMSR("msr_added_turbulence") === FLORIDyn.AddedTurbulence
    @test FLORIDyn.toMSR("msr_eff_wind_speed") === FLORIDyn.EffWind

    # Unknown should throw a helpful error and contain message text
    err = try
        FLORIDyn.toMSR("unknown_kind")
        nothing
    catch e
        e
    end
    @test err !== nothing
    @test isa(err, ErrorException)
    msg = sprint(showerror, err)
    @test occursin("Unknown measurement type", msg)
end

nothing
