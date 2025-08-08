# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test

@testset verbose=true "settings" begin
    @testset "getTurbineData" begin
        # Test single known turbine types
        @testset "Single turbines" begin
            for (name, expected_pos, expected_D) in [
                ("DTU 10MW", [0.0, 0.0, 119.0], 178.4),
                ("DTU 5MW", [0.0, 0.0, 119.0], 178.4),
                ("Senvion 6.2M", [0.0, 0.0, 123.0], 126.0),
                ("V116", [0.0, 0.0, 84.0], 116.0),
                ("V117", [0.0, 0.0, 84.0], 117.0),
                ("V162", [0.0, 0.0, 119.0], 162.0),
                ("GE Haliade X", [0.0, 0.0, 150.0], 220.0)
            ]
                data = getTurbineData([name])
                @test size(data.NacPos) == (1,3)
                @test data.NacPos[1, :] ≈ expected_pos
                @test data.D[1] ≈ expected_D
            end
        end

        # Test multiple turbines
        @testset "Multiple turbines" begin
            names = ["DTU 10MW", "V117", "GE Haliade X"]
            data = getTurbineData(names)
            @test size(data.NacPos) == (3,3)
            @test data.NacPos[1, :] ≈ [0.0, 0.0, 119.0]
            @test data.D[1] ≈ 178.4
            @test data.NacPos[2, :] ≈ [0.0, 0.0, 84.0]
            @test data.D[2] ≈ 117.0
            @test data.NacPos[3, :] ≈ [0.0, 0.0, 150.0]
            @test data.D[3] ≈ 220.0
        end

        # # Test unknown turbine raises error
        # @testset "Unknown turbine error" begin
        #     err = @test_throws ErrorException getTurbineData(["Unknown Turbine"])
        #     println(err.msg)
        #     @test occursin("not known or misspelled", err.msg)
        # end

        # Test empty input returns empty arrays
        @testset "Empty input" begin
            data = getTurbineData(String[])
            @test size(data.NacPos) == (0, 3)
            @test length(data.D) == 0
        end
    end

    @testset "Vis constructor from YAML" begin
        # Test loading Vis from YAML file
        @testset "Load from vis_default.yaml" begin
            vis_file = joinpath(dirname(pathof(FLORIDyn)), "..", "data", "vis_default.yaml")
            vis = Vis(vis_file)
            
            # Test that all fields are loaded correctly
            @test vis.online == false
            @test vis.save == true
            @test vis.print_filenames == false
            @test vis.video_folder == "video"
            @test vis.output_folder == "out"
            @test vis.v_min ≈ 2.0
            @test vis.v_max ≈ 10.0
            @test vis.rel_v_min ≈ 20.0
            @test vis.rel_v_max ≈ 100.0
            @test vis.turb_max ≈ 35.0
            @test vis.up_int == 1
            @test vis.unit_test == false
        end

        # Test computed properties
        @testset "Computed properties" begin
            vis_file = joinpath(dirname(pathof(FLORIDyn)), "..", "data", "vis_default.yaml")
            vis = Vis(vis_file)
            
            # Test video_path property
            video_path = vis.video_path
            @test isa(video_path, String)
            @test endswith(video_path, "video")
            @test isdir(video_path)  # Should be created automatically
            
            # Test output_path property  
            output_path = vis.output_path
            @test isa(output_path, String)
            @test endswith(output_path, "out")
            @test isdir(output_path)  # Should be created automatically
        end

        # Test Base.getproperty method comprehensively
        @testset "Base.getproperty method" begin
            # Create a test Vis struct with custom folder names
            vis = Vis(online=true, save=false, video_folder="test_video", output_folder="test_out")
            
            # Test regular field access
            @testset "Regular field access" begin
                @test vis.online == true
                @test vis.save == false
                @test vis.print_filenames == false
                @test vis.video_folder == "test_video"
                @test vis.output_folder == "test_out"
                @test vis.v_min ≈ 2.0
                @test vis.v_max ≈ 10.0
                @test vis.rel_v_min ≈ 20.0
                @test vis.rel_v_max ≈ 100.0
                @test vis.turb_max ≈ 35.0
                @test vis.up_int == 1
                @test vis.unit_test == false
            end
            
            # Test computed video_path property
            @testset "Computed video_path property" begin
                video_path = vis.video_path
                @test isa(video_path, String)
                @test endswith(video_path, "test_video")
                @test isdir(video_path)  # Directory should be created automatically
                
                # Test path construction based on environment
                if FLORIDyn.isdelftblue()
                    expected_path = joinpath(homedir(), "scratch", "test_video")
                else
                    expected_path = joinpath(pwd(), "test_video")
                end
                @test video_path == expected_path
                
                # Test that accessing multiple times returns the same path
                video_path2 = vis.video_path
                @test video_path == video_path2
            end
            
            # Test computed output_path property
            @testset "Computed output_path property" begin
                output_path = vis.output_path
                @test isa(output_path, String)
                @test endswith(output_path, "test_out")
                @test isdir(output_path)  # Directory should be created automatically
                
                # Test path construction based on environment
                if FLORIDyn.isdelftblue()
                    expected_path = joinpath(homedir(), "scratch", "test_out")
                else
                    expected_path = joinpath(pwd(), "test_out")
                end
                @test output_path == expected_path
                
                # Test that accessing multiple times returns the same path
                output_path2 = vis.output_path
                @test output_path == output_path2
            end
            
            # Test that mkpath creates directories if they don't exist
            @testset "Directory creation" begin
                # Create a unique folder name to ensure it doesn't exist
                unique_folder = "test_unique_$(rand(10000:99999))"
                vis_temp = Vis(online=false, video_folder=unique_folder, output_folder=unique_folder * "_out")
                
                # Before accessing, directories shouldn't exist (in most cases)
                expected_video_path = FLORIDyn.isdelftblue() ? joinpath(homedir(), "scratch", unique_folder) : joinpath(pwd(), unique_folder)
                expected_output_path = FLORIDyn.isdelftblue() ? joinpath(homedir(), "scratch", unique_folder * "_out") : joinpath(pwd(), unique_folder * "_out")
                
                # Access the computed properties - this should create the directories
                video_path = vis_temp.video_path
                output_path = vis_temp.output_path
                
                # Verify directories were created
                @test isdir(video_path)
                @test isdir(output_path)
                @test video_path == expected_video_path
                @test output_path == expected_output_path
                
                # Clean up test directories
                try
                    rm(video_path, recursive=true)
                    rm(output_path, recursive=true)
                catch
                    # Ignore cleanup errors in tests
                end
            end
            
            # Test that getfield still works for accessing underlying fields
            @testset "Direct field access via getfield" begin
                @test getfield(vis, :online) == true
                @test getfield(vis, :save) == false
                @test getfield(vis, :video_folder) == "test_video"
                @test getfield(vis, :output_folder) == "test_out"
                @test getfield(vis, :v_min) ≈ 2.0
            end
            
            # Clean up test directories
            try
                rm(vis.video_path, recursive=true)
                rm(vis.output_path, recursive=true)
            catch
                # Ignore cleanup errors
            end
        end

        # Test error handling for non-existent file
        @testset "Non-existent file error" begin
            @test_throws SystemError Vis("nonexistent_file.yaml")
        end
    end
end
nothing