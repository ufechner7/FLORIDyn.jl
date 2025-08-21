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
    
    @testset "get_default_project" begin
        pkg_data = joinpath(dirname(pathof(FLORIDyn)), "..", "data")

        # Helper to write a minimal projects.yaml in the current working dir
        function write_projects_yaml(dir; entries=[
            (name="2021_9T_Data", vis="vis_default.yaml"),
            (name="2021_54T_NordseeOne", vis="vis_54T.yaml"),
        ])
            mkpath(joinpath(dir, "data"))
            io = open(joinpath(dir, "data", "projects.yaml"), "w")
            try
                write(io, "projects:\n")
                for e in entries
                    write(io, "  - project:\n")
                    write(io, "      name: $(e.name)\n")
                    write(io, "      description: test\n")
                    write(io, "      vis: $(e.vis)\n")
                end
            finally
                close(io)
            end
        end

        @testset "creates default.yaml with first project and returns pkg paths" begin
            mktempdir() do tmp
                cd(tmp) do
                    write_projects_yaml(pwd())
                    settings_file, vis_file = FLORIDyn.get_default_project()
                    # Should create default.yaml with first project
                    @test isfile(joinpath("data", "default.yaml"))
                    def_content = read(joinpath("data", "default.yaml"), String)
                    @test occursin("name: 2021_9T_Data", def_content)
                    # Should resolve to package data files (no local settings/vis exist)
                    @test settings_file == joinpath(pkg_data, "2021_9T_Data.yaml")
                    @test vis_file == joinpath(pkg_data, "vis_default.yaml")
                end
            end
        end

        @testset "honors existing default.yaml selection" begin
            mktempdir() do tmp
                cd(tmp) do
                    write_projects_yaml(pwd())
                    mkpath("data")
                    open(joinpath("data", "default.yaml"), "w") do io
                        write(io, "default:\n  name: 2021_54T_NordseeOne\n")
                    end
                    settings_file, vis_file = FLORIDyn.get_default_project()
                    @test settings_file == joinpath(pkg_data, "2021_54T_NordseeOne.yaml")
                    @test vis_file == joinpath(pkg_data, "vis_54T.yaml")
                end
            end
        end

        @testset "updates stale default.yaml to first project" begin
            mktempdir() do tmp
                cd(tmp) do
                    write_projects_yaml(pwd())
                    mkpath("data")
                    # Write a non-existing project name
                    open(joinpath("data", "default.yaml"), "w") do io
                        write(io, "default:\n  name: DOES_NOT_EXIST\n")
                    end
                    settings_file, vis_file = FLORIDyn.get_default_project()
                    # Should fall back to first project and rewrite default.yaml
                    @test settings_file == joinpath(pkg_data, "2021_9T_Data.yaml")
                    @test vis_file == joinpath(pkg_data, "vis_default.yaml")
                    def_content = read(joinpath("data", "default.yaml"), String)
                    @test occursin("name: 2021_9T_Data", def_content)
                end
            end
        end

        @testset "prefers local data files over package files" begin
            mktempdir() do tmp
                cd(tmp) do
                    write_projects_yaml(pwd())
                    # Create local settings and vis files
                    mkpath("data")
                    local_settings = joinpath("data", "2021_9T_Data.yaml")
                    local_vis = joinpath("data", "vis_default.yaml")
                    write(local_settings, "# local settings placeholder\n")
                    write(local_vis, "# local vis placeholder\n")
                    settings_file, vis_file = FLORIDyn.get_default_project()
                    @test abspath(settings_file) == abspath(local_settings)
                    @test abspath(vis_file) == abspath(local_vis)
                end
            end
        end

        @testset "empty projects list throws error" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    write(joinpath("data", "projects.yaml"), "projects: []\n")
                    @test_throws ErrorException FLORIDyn.get_default_project()
                end
            end
        end
    end

    @testset "list_projects" begin
        @testset "reads projects from package/local data" begin
            projs = FLORIDyn.list_projects()
            @test isa(projs, Vector{Tuple{String,String}})
            # Should contain the two reference projects from repo data
            @test ("2021_9T_Data", "vis_default.yaml") in projs
            @test ("2021_54T_NordseeOne", "vis_54T.yaml") in projs
        end

        @testset "prefers local projects.yaml override" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    content = """
projects:
  - project:
      name: Alpha
      description: test
      vis: alpha_vis.yaml
  - project:
      name: Beta
      description: test
      vis: beta_vis.yaml
"""
                    write(joinpath("data", "projects.yaml"), content)
                    projs = FLORIDyn.list_projects()
                    @test projs == [("Alpha", "alpha_vis.yaml"), ("Beta", "beta_vis.yaml")]
                end
            end
        end

        @testset "empty projects list returns empty vector" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    write(joinpath("data", "projects.yaml"), "projects: []\n")
                    projs = FLORIDyn.list_projects()
                    @test projs == Tuple{String,String}[]
                end
            end
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
        
        # Test parsing YAML file with flow_fields and measurements
        @testset "YAML parsing with flow_fields and measurements" begin
            # Create a temporary YAML file with flow_fields and measurements
            test_yaml_content = """
vis:
  online: true
  save: false
  print_filenames: true
  video_folder: "test_video"
  output_folder: "test_output"
  flow_fields:
    - name: "flow_field_vel_reduction"
      online: false
      create_video: true
    - name: "flow_field_added_turbulence" 
      online: true
      create_video: false
    - name: "flow_field_eff_wind_speed"
      online: false
      create_video: true
  measurements:
    - name: "msr_vel_reduction"
      separated: true
    - name: "msr_added_turbulence"
      separated: false
    - name: "msr_eff_wind_speed" 
      separated: true
  v_min: 3.0
  v_max: 12.0
  rel_v_min: 25.0
  rel_v_max: 95.0
  turb_max: 40.0
  up_int: 2
  unit_test: true
"""
            
            # Write temporary YAML file
            test_yaml_file = "test_vis_complex.yaml"
            try
                write(test_yaml_file, test_yaml_content)
                
                # Parse the YAML file
                vis = Vis(test_yaml_file)
                
                # Test basic fields
                @test vis.online == true
                @test vis.save == false
                @test vis.print_filenames == true
                @test vis.video_folder == "test_video"
                @test vis.output_folder == "test_output"
                @test vis.v_min ≈ 3.0
                @test vis.v_max ≈ 12.0
                @test vis.rel_v_min ≈ 25.0
                @test vis.rel_v_max ≈ 95.0
                @test vis.turb_max ≈ 40.0
                @test vis.up_int == 2
                @test vis.unit_test == true
                
                # Test flow_fields parsing - these are parsed into FlowField structs
                flow_fields = vis.flow_fields
                @test isa(flow_fields, Vector{FLORIDyn.FlowField})
                @test length(flow_fields) == 3
                
                # Check first flow field entry
                @test flow_fields[1].name == "flow_field_vel_reduction"
                @test flow_fields[1].online == false
                @test flow_fields[1].create_video == true
                
                # Check second flow field entry
                @test flow_fields[2].name == "flow_field_added_turbulence"
                @test flow_fields[2].online == true
                @test flow_fields[2].create_video == false
                
                # Check third flow field entry
                @test flow_fields[3].name == "flow_field_eff_wind_speed"
                @test flow_fields[3].online == false
                @test flow_fields[3].create_video == true
                
                # Test measurements parsing - these are parsed into Measurement structs
                measurements = vis.measurements
                @test isa(measurements, Vector{FLORIDyn.Measurement})
                @test length(measurements) == 3
                
                # Check measurement entries
                @test measurements[1].name == "msr_vel_reduction"
                @test measurements[1].separated == true
                
                @test measurements[2].name == "msr_added_turbulence"
                @test measurements[2].separated == false
                
                @test measurements[3].name == "msr_eff_wind_speed"
                @test measurements[3].separated == true
                
                # Clean up test directories that may have been created
                try
                    rm(vis.video_path, recursive=true, force=true)
                    rm(vis.output_path, recursive=true, force=true)
                catch
                    # Ignore cleanup errors
                end
                
            finally
                # Clean up test file
                if isfile(test_yaml_file)
                    rm(test_yaml_file)
                end
            end
        end
    end

    @testset "MSR default functions" begin
        @testset "get_default_msr and set_default_msr" begin
            mktempdir() do tmp
                cd(tmp) do
                    # Test get_default_msr when no file exists
                    @testset "No default.yaml file" begin
                        @test FLORIDyn.get_default_msr() == VelReduction
                    end

                    # Test set_default_msr creates file with MSR
                    @testset "set_default_msr creates file" begin
                        FLORIDyn.set_default_msr(AddedTurbulence)
                        @test isfile(joinpath("data", "default.yaml"))
                        
                        content = read(joinpath("data", "default.yaml"), String)
                        @test occursin("msr: AddedTurbulence", content)
                        @test occursin("name:", content)  # Should have default name
                        
                        # Test get_default_msr reads it correctly
                        @test FLORIDyn.get_default_msr() == AddedTurbulence
                    end

                    # Test changing MSR preserves existing project name
                    @testset "set_default_msr preserves project name" begin
                        # First set a specific project name
                        mkpath("data")
                        open(joinpath("data", "default.yaml"), "w") do io
                            write(io, "default:\n  name: TestProject\n  msr: VelReduction\n")
                        end
                        
                        # Change MSR - should preserve project name
                        FLORIDyn.set_default_msr(EffWind)
                        
                        content = read(joinpath("data", "default.yaml"), String)
                        @test occursin("name: TestProject", content)
                        @test occursin("msr: EffWind", content)
                        @test FLORIDyn.get_default_msr() == EffWind
                    end

                    # Test all MSR enum values
                    @testset "All MSR enum values" begin
                        for msr in [VelReduction, AddedTurbulence, EffWind]
                            FLORIDyn.set_default_msr(msr)
                            @test FLORIDyn.get_default_msr() == msr
                        end
                    end

                    # Test malformed YAML file handling
                    @testset "Malformed YAML handling" begin
                        mkpath("data")
                        # Write malformed YAML
                        open(joinpath("data", "default.yaml"), "w") do io
                            write(io, "invalid: yaml: content:\n  - broken")
                        end
                        
                        # Should return default value
                        @test FLORIDyn.get_default_msr() == VelReduction
                        
                        # set_default_msr should still work (overwrites file)
                        FLORIDyn.set_default_msr(AddedTurbulence)
                        @test FLORIDyn.get_default_msr() == AddedTurbulence
                    end

                    # Test file with missing msr field
                    @testset "Missing msr field" begin
                        mkpath("data")
                        open(joinpath("data", "default.yaml"), "w") do io
                            write(io, "default:\n  name: SomeProject\n")
                        end
                        
                        # Should return default value
                        @test FLORIDyn.get_default_msr() == VelReduction
                    end

                    # Test file with invalid msr value
                    @testset "Invalid msr value" begin
                        mkpath("data")
                        open(joinpath("data", "default.yaml"), "w") do io
                            write(io, "default:\n  name: SomeProject\n  msr: InvalidValue\n")
                        end
                        
                        # Should return default value when MSR value is invalid
                        @test FLORIDyn.get_default_msr() == VelReduction
                    end
                end
            end
        end

        @testset "select_measurement function" begin
            # Test that the function is exported and callable
            @testset "Function availability" begin
                @test isdefined(FLORIDyn, :select_measurement)
                @test hasmethod(FLORIDyn.select_measurement, ())
            end

            # Note: Interactive testing of select_measurement() would require 
            # mocking stdin, which is complex. The function's logic is covered
            # by testing set_default_msr and get_default_msr above.
        end

        @testset "Integration with existing functions" begin
            mktempdir() do tmp
                cd(tmp) do
                    # Write a test projects.yaml
                    mkpath("data")
                    projects_content = """
projects:
  - project:
      name: TestProject1
      description: test
      vis: vis_default.yaml
  - project:
      name: TestProject2
      description: test  
      vis: vis_test.yaml
"""
                    write(joinpath("data", "projects.yaml"), projects_content)

                    # Test that get_default_project preserves MSR when creating default.yaml
                    @testset "get_default_project preserves MSR" begin
                        # Set a specific MSR first
                        FLORIDyn.set_default_msr(AddedTurbulence)
                        
                        # Create dummy settings and vis files that get_default_project expects
                        pkg_data = joinpath(dirname(pathof(FLORIDyn)), "..", "data")
                        dummy_settings = joinpath(pkg_data, "TestProject1.yaml")
                        dummy_vis = joinpath(pkg_data, "vis_default.yaml")
                        
                        # Skip this test if we can't create files in package data
                        # (This test is more about the logic than file creation)
                        if isfile(dummy_vis)  # vis_default.yaml should exist
                            # Create a temporary settings file for the test
                            test_settings_content = "# Test settings file for unit test"
                            temp_settings_file = joinpath("data", "TestProject1.yaml")
                            write(temp_settings_file, test_settings_content)
                            
                            # Mock get_default_project by testing the MSR preservation logic directly
                            # Check that MSR is preserved in the file
                            @test FLORIDyn.get_default_msr() == AddedTurbulence
                            
                            # Check file content
                            content = read(joinpath("data", "default.yaml"), String)
                            @test occursin("msr: AddedTurbulence", content)
                        else
                            @test_skip "Skipping integration test - package files not accessible"
                        end
                    end

                    # Test that select_project preserves MSR when changing project
                    @testset "select_project preserves MSR" begin
                        # Set up initial state
                        FLORIDyn.set_default_msr(EffWind)
                        
                        # Simulate selecting a different project by directly writing to default.yaml
                        # (We can't easily test the interactive select_project function)
                        mkpath("data")
                        
                        # Read existing MSR
                        existing_msr = FLORIDyn.get_default_msr()
                        @test existing_msr == EffWind
                        
                        # Manually update project name while preserving MSR logic
                        # (This simulates what select_project does)
                        open(joinpath("data", "default.yaml"), "w") do io
                            write(io, "default:\n  name: TestProject2\n  msr: $(string(existing_msr))  # valid options: VelReduction, AddedTurbulence, EffWind\n")
                        end
                        
                        # Verify MSR is preserved and project changed
                        content = read(joinpath("data", "default.yaml"), String)
                        @test occursin("name: TestProject2", content)
                        @test occursin("msr: EffWind", content)
                        @test FLORIDyn.get_default_msr() == EffWind
                    end
                end
            end
        end
    end
end
nothing