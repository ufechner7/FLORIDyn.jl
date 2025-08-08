# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Test
using FLORIDyn

@testset "Measurement Struct Tests" begin
    @testset "Measurement Construction" begin
        # Test basic construction
        m1 = Measurement("test_measurement", true)
        @test m1.name == "test_measurement"
        @test m1.separated == true
        
        # Test default separated value
        m2 = Measurement("another_measurement")
        @test m2.name == "another_measurement"
        @test m2.separated == false
        
        # Test with keyword arguments
        m3 = Measurement(name="kw_measurement", separated=true)
        @test m3.name == "kw_measurement"
        @test m3.separated == true
    end
    
    @testset "parse_measurements Function" begin
        # Test simple string format
        simple_yaml = ["measurements_vel_reduction", "measurements_added_turbulence"]
        measurements = parse_measurements(simple_yaml)
        
        @test length(measurements) == 2
        @test measurements[1].name == "measurements_vel_reduction"
        @test measurements[1].separated == false
        @test measurements[2].name == "measurements_added_turbulence"
        @test measurements[2].separated == false
        
        # Test dictionary format
        dict_yaml = [
            Dict("name" => "measurements_vel_reduction", "separated" => true),
            Dict("name" => "measurements_added_turbulence", "separated" => false)
        ]
        measurements = parse_measurements(dict_yaml)
        
        @test length(measurements) == 2
        @test measurements[1].name == "measurements_vel_reduction"
        @test measurements[1].separated == true
        @test measurements[2].name == "measurements_added_turbulence"
        @test measurements[2].separated == false
        
        # Test mixed format
        mixed_yaml = [
            "measurements_vel_reduction",
            Dict("name" => "measurements_added_turbulence", "separated" => true),
            "measurements_eff_wind_speed"
        ]
        measurements = parse_measurements(mixed_yaml)
        
        @test length(measurements) == 3
        @test measurements[1].name == "measurements_vel_reduction"
        @test measurements[1].separated == false
        @test measurements[2].name == "measurements_added_turbulence"
        @test measurements[2].separated == true
        @test measurements[3].name == "measurements_eff_wind_speed"
        @test measurements[3].separated == false
        
        # Test with incomplete dictionary
        incomplete_yaml = [
            Dict("name" => "complete_measurement", "separated" => true),
            Dict("separated" => true),  # Missing name
            Dict("name" => "name_only")  # Missing separated
        ]
        measurements = parse_measurements(incomplete_yaml)
        
        @test length(measurements) == 2  # Only complete entries should be parsed
        @test measurements[1].name == "complete_measurement"
        @test measurements[1].separated == true
        @test measurements[2].name == "name_only"
        @test measurements[2].separated == false  # Default value
        
        # Test empty input
        empty_measurements = parse_measurements([])
        @test length(empty_measurements) == 0
    end
    
    @testset "parse_measurements with Dict input" begin
        # Test parsing from vis configuration dictionary
        vis_config = Dict(
            "measurements" => [
                Dict("name" => "measurements_vel_reduction", "separated" => true),
                "measurements_added_turbulence"
            ]
        )
        
        measurements = parse_measurements(vis_config)
        @test length(measurements) == 2
        @test measurements[1].name == "measurements_vel_reduction"
        @test measurements[1].separated == true
        @test measurements[2].name == "measurements_added_turbulence"
        @test measurements[2].separated == false
        
        # Test with missing measurements key
        empty_config = Dict("other_key" => "value")
        measurements = parse_measurements(empty_config)
        @test length(measurements) == 0
    end
end
