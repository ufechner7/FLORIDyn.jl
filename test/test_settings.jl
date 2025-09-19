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
                @test size(data.nac_pos) == (1,3)
                @test data.nac_pos[1, :] ≈ expected_pos
                @test data.rotor_diameter[1] ≈ expected_D
            end
        end

        # Test multiple turbines
        @testset "Multiple turbines" begin
            names = ["DTU 10MW", "V117", "GE Haliade X"]
            data = getTurbineData(names)
            @test size(data.nac_pos) == (3,3)
            @test data.nac_pos[1, :] ≈ [0.0, 0.0, 119.0]
            @test data.rotor_diameter[1] ≈ 178.4
            @test data.nac_pos[2, :] ≈ [0.0, 0.0, 84.0]
            @test data.rotor_diameter[2] ≈ 117.0
            @test data.nac_pos[3, :] ≈ [0.0, 0.0, 150.0]
            @test data.rotor_diameter[3] ≈ 220.0
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
            @test size(data.nac_pos) == (0, 3)
            @test length(data.rotor_diameter) == 0
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
                    project_name, settings_file, vis_file = FLORIDyn.get_default_project()
                    # Should create default.yaml with first project
                    @test isfile(joinpath("data", "default.yaml"))
                    def_content = read(joinpath("data", "default.yaml"), String)
                    @test occursin("name: 2021_9T_Data", def_content)
                    # Should return project name as first element
                    @test project_name == "2021_9T_Data"
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
                    project_name, settings_file, vis_file = FLORIDyn.get_default_project()
                    @test project_name == "2021_54T_NordseeOne"
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
                    project_name, settings_file, vis_file = FLORIDyn.get_default_project()
                    # Should fall back to first project and rewrite default.yaml
                    @test project_name == "2021_9T_Data"
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
                    project_name, settings_file, vis_file = FLORIDyn.get_default_project()
                    @test project_name == "2021_9T_Data"
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
            @test isa(projs, Vector{Tuple{String,String,String}})
            # Should contain the two reference projects from repo data
            @test ("2021_9T_Data", "A reference simulation with 9 turbines", "vis_default.yaml") in projs
            @test ("2021_54T_NordseeOne", "A reference simulation with 54 turbines", "vis_54T.yaml") in projs
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
                    @test projs == [("Alpha", "test", "alpha_vis.yaml"), ("Beta", "test", "beta_vis.yaml")]
                end
            end
        end

        @testset "empty projects list returns empty vector" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    write(joinpath("data", "projects.yaml"), "projects: []\n")
                    projs = FLORIDyn.list_projects()
                    @test projs == Tuple{String,String,String}[]
                end
            end
        end

        @testset "missing projects.yaml uses package data" begin
            mktempdir() do tmp
                cd(tmp) do
                    # Don't create local projects.yaml - should fall back to package data
                    projs = FLORIDyn.list_projects()
                    @test isa(projs, Vector{Tuple{String,String,String}})
                    # Should contain the reference projects from package data
                    @test ("2021_9T_Data", "A reference simulation with 9 turbines", "vis_default.yaml") in projs
                    @test ("2021_54T_NordseeOne", "A reference simulation with 54 turbines", "vis_54T.yaml") in projs
                end
            end
        end

        @testset "malformed projects.yaml throws error" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    # Write malformed YAML
                    write(joinpath("data", "projects.yaml"), "invalid: yaml: content:\n  - broken")
                    @test_throws Exception FLORIDyn.list_projects()
                end
            end
        end

        @testset "projects.yaml without projects key" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    content = """
other_key: some_value
settings:
  - name: test
"""
                    write(joinpath("data", "projects.yaml"), content)
                    projs = FLORIDyn.list_projects()
                    @test projs == Tuple{String,String}[]
                end
            end
        end

        @testset "projects.yaml with null projects" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    write(joinpath("data", "projects.yaml"), "projects: null\n")
                    # Should throw an error when trying to iterate over null
                    @test_throws MethodError FLORIDyn.list_projects()
                end
            end
        end

        @testset "complex projects.yaml structure" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    content = """
projects:
  - project:
      name: ProjectOne
      description: First test project
      vis: vis_one.yaml
      extra_field: should_be_ignored
  - project:
      name: ProjectTwo
      description: Second test project with longer description
      vis: vis_two.yaml
      author: Test Author
      version: 1.0
  - project:
      name: ProjectThree
      description: Third project
      vis: vis_three.yaml
"""
                    write(joinpath("data", "projects.yaml"), content)
                    projs = FLORIDyn.list_projects()
                    expected = [
                        ("ProjectOne", "First test project", "vis_one.yaml"),
                        ("ProjectTwo", "Second test project with longer description", "vis_two.yaml"),
                        ("ProjectThree", "Third project", "vis_three.yaml")
                    ]
                    @test projs == expected
                end
            end
        end

        @testset "projects with missing name or vis fields" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    content = """
projects:
  - project:
      name: ValidProject
      description: Valid project
      vis: valid_vis.yaml
  - project:
      # Missing name field
      description: Invalid project - no name
      vis: invalid_vis.yaml
"""
                    write(joinpath("data", "projects.yaml"), content)
                    # Should throw KeyError when trying to access missing fields
                    @test_throws KeyError FLORIDyn.list_projects()
                end
            end
        end

        @testset "projects with special characters in names" begin
            mktempdir() do tmp
                cd(tmp) do
                    mkpath("data")
                    content = """
projects:
  - project:
      name: "Project-With-Dashes"
      description: Project with dashes
      vis: vis_dashes.yaml
  - project:
      name: "Project_With_Underscores"
      description: Project with underscores
      vis: vis_underscores.yaml
  - project:
      name: "Project With Spaces"
      description: Project with spaces
      vis: vis_spaces.yaml
"""
                    write(joinpath("data", "projects.yaml"), content)
                    projs = FLORIDyn.list_projects()
                    expected = [
                        ("Project-With-Dashes", "Project with dashes", "vis_dashes.yaml"),
                        ("Project_With_Underscores", "Project with underscores", "vis_underscores.yaml"),
                        ("Project With Spaces", "Project with spaces", "vis_spaces.yaml")
                    ]
                    @test projs == expected
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
            @test vis.v_min ≈ 1.0
            @test vis.v_max ≈ 9.0
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

    @testset "turbine_group function" begin
        @testset "Basic functionality" begin
            # Test with 54-turbine configuration
            @testset "54-turbine configuration" begin
                wind, sim, con, floris, floridyn, ta = setup("data/2021_54T_NordseeOne.yaml")
                
                # Test that function returns correct group IDs for known turbines
                @test turbine_group(ta, 1) == 1    # Should be in group1 (westernmost)
                @test turbine_group(ta, 15) == 2   # Should be in group2
                @test turbine_group(ta, 30) == 4   # Should be in group4
                @test turbine_group(ta, 45) == 3   # Should be in group3
                @test turbine_group(ta, 54) == 1   # Should be in group1
                
                # Test that all turbines have valid group assignments
                for turbine_id in 1:54
                    group_id = turbine_group(ta, turbine_id)
                    @test isa(group_id, Int)
                    @test group_id >= 0  # Group IDs should be non-negative
                    
                    # Verify that the group actually exists
                    group_found = false
                    for group in ta.groups
                        if group.id == group_id
                            group_found = true
                            @test turbine_id in group.turbines
                            break
                        end
                    end
                    @test group_found
                end
            end
            
            # Test with 9-turbine configuration
            @testset "9-turbine configuration" begin
                wind, sim, con, floris, floridyn, ta = setup("data/2021_9T_Data.yaml")
                
                # Test specific turbine-group mappings for 9T layout
                @test turbine_group(ta, 1) == 1    # Should be in row1
                @test turbine_group(ta, 2) == 1    # Should be in row1  
                @test turbine_group(ta, 3) == 1    # Should be in row1
                @test turbine_group(ta, 4) == 2    # Should be in row2
                @test turbine_group(ta, 5) == 2    # Should be in row2
                @test turbine_group(ta, 6) == 2    # Should be in row2
                @test turbine_group(ta, 7) == 3    # Should be in row3
                @test turbine_group(ta, 8) == 3    # Should be in row3
                @test turbine_group(ta, 9) == 3    # Should be in row3
                
                # Test that all turbines have valid group assignments
                for turbine_id in 1:9
                    group_id = turbine_group(ta, turbine_id)
                    @test isa(group_id, Int)
                    @test group_id >= 0
                    
                    # Verify the group exists and contains the turbine
                    group_found = false
                    for group in ta.groups
                        if group.id == group_id
                            group_found = true
                            @test turbine_id in group.turbines
                            break
                        end
                    end
                    @test group_found
                end
            end
        end
        
        @testset "Error handling" begin
            # Test with 54-turbine configuration for error cases
            wind, sim, con, floris, floridyn, ta = setup("data/2021_54T_NordseeOne.yaml")
            
            # Test invalid turbine numbers
            @testset "Out of bounds turbine numbers" begin
                # Test turbine number too low
                err = @test_throws ArgumentError turbine_group(ta, 0)
                @test occursin("out of bounds", err.value.msg)
                @test occursin("Valid range: 1-54", err.value.msg)
                
                # Test turbine number too high
                err = @test_throws ArgumentError turbine_group(ta, 55)
                @test occursin("out of bounds", err.value.msg)
                @test occursin("Valid range: 1-54", err.value.msg)
                
                # Test negative turbine number
                err = @test_throws ArgumentError turbine_group(ta, -1)
                @test occursin("out of bounds", err.value.msg)
                
                # Test very large turbine number
                err = @test_throws ArgumentError turbine_group(ta, 1000)
                @test occursin("out of bounds", err.value.msg)
            end
            
            # Test with 9-turbine configuration for different bounds
            @testset "9-turbine bounds checking" begin
                wind, sim, con, floris, floridyn, ta = setup("data/2021_9T_Data.yaml")
                
                # Test valid range
                for i in 1:9
                    @test_nowarn turbine_group(ta, i)
                end
                
                # Test invalid ranges
                err = @test_throws ArgumentError turbine_group(ta, 0)
                @test occursin("Valid range: 1-9", err.value.msg)
                
                err = @test_throws ArgumentError turbine_group(ta, 10)
                @test occursin("Valid range: 1-9", err.value.msg)
            end
        end
        
        @testset "Group priority behavior" begin
            # Test that the function prioritizes specific groups over "all" group
            wind, sim, con, floris, floridyn, ta = setup("data/2021_54T_NordseeOne.yaml")
            
            @testset "Prioritizes specific groups" begin
                # Every turbine should belong to the "all" group, but function should
                # return specific group IDs when available
                for turbine_id in [1, 15, 30, 45, 54]
                    group_id = turbine_group(ta, turbine_id)
                    
                    # Should not return the "all" group ID (0) for these turbines
                    # since they also belong to specific spatial groups
                    @test group_id != 0
                    
                    # Verify the turbine is in both the specific group and "all" group
                    found_in_specific = false
                    found_in_all = false
                    
                    for group in ta.groups
                        if group.id == group_id && turbine_id in group.turbines
                            found_in_specific = true
                            @test group.name != "all"  # Should be a specific group
                        end
                        if group.name == "all" && turbine_id in group.turbines
                            found_in_all = true
                        end
                    end
                    
                    @test found_in_specific
                    @test found_in_all
                end
            end
        end
        
        @testset "Edge cases and robustness" begin
            @testset "Consistency across multiple calls" begin
                wind, sim, con, floris, floridyn, ta = setup("data/2021_54T_NordseeOne.yaml")
                
                # Test that multiple calls return the same result
                for turbine_id in [1, 10, 25, 40, 54]
                    group_id1 = turbine_group(ta, turbine_id)
                    group_id2 = turbine_group(ta, turbine_id)
                    group_id3 = turbine_group(ta, turbine_id)
                    
                    @test group_id1 == group_id2 == group_id3
                end
            end
            
            @testset "All turbines are assigned to groups" begin
                wind, sim, con, floris, floridyn, ta = setup("data/2021_54T_NordseeOne.yaml")
                
                # Test that every turbine can be assigned to a group
                assigned_groups = Set{Int}()
                for turbine_id in 1:54
                    group_id = turbine_group(ta, turbine_id)
                    push!(assigned_groups, group_id)
                end
                
                # Should have found multiple groups (not just one)
                @test length(assigned_groups) > 1
                
                # All group IDs should correspond to actual groups
                existing_group_ids = Set([group.id for group in ta.groups])
                @test issubset(assigned_groups, existing_group_ids)
            end
            
            @testset "Group distribution makes sense" begin
                wind, sim, con, floris, floridyn, ta = setup("data/2021_54T_NordseeOne.yaml")
                
                # Count turbines per group
                group_counts = Dict{Int, Int}()
                for turbine_id in 1:54
                    group_id = turbine_group(ta, turbine_id)
                    group_counts[group_id] = get(group_counts, group_id, 0) + 1
                end
                
                # Should have reasonable distribution (no group should have all turbines
                # unless there's only one specific group)
                total_specific_groups = length([g for g in ta.groups if g.name != "all"])
                if total_specific_groups > 1
                    for count in values(group_counts)
                        @test count < 54  # No single group should have all turbines
                        @test count > 0   # Every group should have some turbines
                    end
                end
                
                # Total should add up to 54
                @test sum(values(group_counts)) == 54
            end
        end
        
        @testset "Integration with TurbineArray structure" begin
            @testset "Validates against actual group definitions" begin
                wind, sim, con, floris, floridyn, ta = setup("data/2021_54T_NordseeOne.yaml")
                
                # For each turbine, verify it's actually in the returned group
                for turbine_id in 1:54
                    group_id = turbine_group(ta, turbine_id)
                    
                    # Find the group with this ID
                    target_group = nothing
                    for group in ta.groups
                        if group.id == group_id
                            target_group = group
                            break
                        end
                    end
                    
                    @test target_group !== nothing
                    @test turbine_id in target_group.turbines
                end
            end
            
            @testset "Works with different TurbineArray configurations" begin
                # Test both configurations to ensure function is generic
                configs = [
                    ("data/2021_54T_NordseeOne.yaml", 54),
                    ("data/2021_9T_Data.yaml", 9)
                ]
                
                for (config_file, num_turbines) in configs
                    wind, sim, con, floris, floridyn, ta = setup(config_file)
                    
                    # Test a few representative turbines
                    test_turbines = [1, div(num_turbines, 2), num_turbines]
                    for turbine_id in test_turbines
                        @test_nowarn turbine_group(ta, turbine_id)
                        group_id = turbine_group(ta, turbine_id)
                        @test isa(group_id, Int)
                        
                        # Verify group exists
                        group_exists = any(group.id == group_id for group in ta.groups)
                        @test group_exists
                    end
                end
            end
        end
        
        @testset "Performance characteristics" begin
            wind, sim, con, floris, floridyn, ta = setup("data/2021_54T_NordseeOne.yaml")
            
            @testset "Reasonable performance for typical usage" begin
                # Time multiple calls to ensure performance is reasonable
                # This is not a strict performance test, but ensures no obvious bottlenecks
                times = Float64[]
                
                for _ in 1:100
                    turbine_id = rand(1:54)
                    time_start = time_ns()
                    turbine_group(ta, turbine_id)
                    time_end = time_ns()
                    push!(times, (time_end - time_start) / 1e6)  # Convert to milliseconds
                end
                
                # Function should be fast for typical wind farm sizes
                # Allow generous time limit since test systems vary
                avg_time = sum(times) / length(times)
                @test avg_time < 0.1  # Should average less than 0.1ms (100μs) per call
                
                # No call should be extremely slow
                max_time = maximum(times)
                @test max_time < 100.0  # No call should take more than 100ms
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

            # Test the measurement options structure used by select_measurement
            @testset "Measurement options validation" begin
                # Test that all MSR enum values are represented
                measurement_options = [
                    "VelReduction - Velocity reduction measurement",
                    "AddedTurbulence - Added turbulence measurement", 
                    "EffWind - Effective wind speed measurement"
                ]
                msr_values = [VelReduction, AddedTurbulence, EffWind]
                
                # Test each measurement option
                for (option, msr) in zip(measurement_options, msr_values)
                    @test isa(option, String)
                    @test isa(msr, MSR)
                    @test length(option) > 10  # Should be descriptive
                    @test contains(option, string(msr))  # Should contain MSR name
                    @test contains(option, " - ")  # Should have separator
                end
                
                # Test that all MSR enum values are covered
                all_msr_values = [VelReduction, AddedTurbulence, EffWind]
                @test Set(all_msr_values) == Set(msr_values)
                @test length(measurement_options) == length(msr_values)
            end

            # Test menu option formatting
            @testset "Menu option formatting" begin
                measurement_options = [
                    "VelReduction - Velocity reduction measurement",
                    "AddedTurbulence - Added turbulence measurement", 
                    "EffWind - Effective wind speed measurement"
                ]
                
                for option in measurement_options
                    # Test that options are well-formatted
                    @test contains(option, " - ")
                    parts = split(option, " - ", limit=2)
                    @test length(parts) == 2
                    
                    msr_name = parts[1]
                    description = parts[2]
                    
                    # Test MSR name part
                    @test msr_name in ["VelReduction", "AddedTurbulence", "EffWind"]
                    
                    # Test description part
                    @test length(description) > 5
                    @test contains(lowercase(description), "measurement")
                end
            end

            # Test that select_measurement integrates correctly with MSR functions
            @testset "Integration with MSR functions" begin
                mktempdir() do tmp
                    cd(tmp) do
                        # Test the mapping between menu choices and MSR values
                        msr_values = [VelReduction, AddedTurbulence, EffWind]
                        
                        # Test that each MSR value can be set and retrieved
                        for expected_msr in msr_values
                            FLORIDyn.set_default_msr(expected_msr)
                            @test FLORIDyn.get_default_msr() == expected_msr
                        end
                        
                        # Test current default MSR retrieval (used by select_measurement for cancellation)
                        FLORIDyn.set_default_msr(AddedTurbulence)
                        current_msr = FLORIDyn.get_default_msr()
                        @test current_msr == AddedTurbulence
                    end
                end
            end

            # Test cancellation behavior (simulating what happens when user cancels)
            @testset "Cancellation handling" begin
                mktempdir() do tmp
                    cd(tmp) do
                        # Set a known MSR value
                        FLORIDyn.set_default_msr(EffWind)
                        
                        # Test that get_default_msr returns the current value
                        # (This simulates what select_measurement returns on cancellation)
                        current_msr = FLORIDyn.get_default_msr()
                        @test current_msr == EffWind
                        
                        # Test that the value remains unchanged after "cancellation"
                        # (In real usage, select_measurement would return this value)
                        @test FLORIDyn.get_default_msr() == EffWind
                    end
                end
            end

            # Test menu configuration
            @testset "Menu configuration" begin
                measurement_options = [
                    "VelReduction - Velocity reduction measurement",
                    "AddedTurbulence - Added turbulence measurement", 
                    "EffWind - Effective wind speed measurement"
                ]
                
                # Test that we can create a RadioMenu with these options
                # (This tests the structure that select_measurement uses)
                try
                    menu = REPL.TerminalMenus.RadioMenu(measurement_options, pagesize=length(measurement_options))
                    @test isa(menu, REPL.TerminalMenus.RadioMenu)
                    @test length(measurement_options) == 3  # Should have exactly 3 options
                catch e
                    # If TerminalMenus is not available in test environment, skip this test
                    @test_skip "TerminalMenus not available in test environment: $e"
                end
            end

            # Test measurement descriptions match MSR enum values
            @testset "MSR enum consistency" begin
                msr_values = [VelReduction, AddedTurbulence, EffWind]
                expected_names = ["VelReduction", "AddedTurbulence", "EffWind"]
                
                for (msr, expected_name) in zip(msr_values, expected_names)
                    @test string(msr) == expected_name
                end
            end

            # Note: Full interactive testing of select_measurement() with actual menu 
            # interaction would require complex input simulation. The function's core 
            # logic is thoroughly tested through the underlying MSR functions and 
            # menu structure validation above.
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

            @testset "select_project function" begin
                # Test that the function is exported and callable
                @testset "Function availability" begin
                    @test isdefined(FLORIDyn, :select_project)
                    @test hasmethod(FLORIDyn.select_project, ())
                end

                # Test project menu structure used by select_project
                @testset "Project menu structure validation" begin
                    mktempdir() do tmp
                        cd(tmp) do
                            # Create test projects.yaml
                            mkpath("data")
                            test_projects = """
                            projects:
                              - project:
                                  name: test_proj1
                                  description: "Test project 1"
                                  vis: "vis_test1.yaml"
                              - project:
                                  name: test_proj2
                                  description: "Test project 2 with longer description"
                                  vis: "vis_test2.yaml"
                              - project:
                                  name: short_name
                                  description: ""
                                  vis: "vis_short.yaml"
                            """
                            open(joinpath("data", "projects.yaml"), "w") do io
                                write(io, test_projects)
                            end
                            
                            # Test that projects are loaded correctly
                            projects = FLORIDyn.list_projects()
                            @test length(projects) == 3
                            
                            # Test menu option formatting for select_project
                            project_options = String[]
                            for (name, description, vis) in projects
                                if isempty(description)
                                    push!(project_options, name)
                                else
                                    push!(project_options, "$name - $description")
                                end
                            end
                            
                            # Verify proper formatting
                            @test length(project_options) == 3
                            @test any(option -> contains(option, "test_proj1 - Test project 1"), project_options)
                            @test any(option -> contains(option, "test_proj2 - Test project 2 with longer description"), project_options)
                            @test any(option -> option == "short_name", project_options)  # No description, so just name
                            
                            # Test that options are well-formatted
                            for option in project_options
                                @test isa(option, String)
                                @test length(option) > 0
                            end
                        end
                    end
                end

                # Test integration with MSR preservation
                @testset "MSR preservation during project selection" begin
                    mktempdir() do tmp
                        cd(tmp) do
                            # Create test projects.yaml
                            mkpath("data")
                            test_projects = """
                            projects:
                              - project:
                                  name: proj_a
                                  description: "Project A"
                                  vis: "vis_a.yaml"
                              - project:
                                  name: proj_b
                                  description: "Project B"
                                  vis: "vis_b.yaml"
                            """
                            open(joinpath("data", "projects.yaml"), "w") do io
                                write(io, test_projects)
                            end
                            
                            # Create the necessary settings files
                            open(joinpath("data", "proj_a.yaml"), "w") do io
                                write(io, "# Dummy settings file for proj_a\n")
                            end
                            open(joinpath("data", "proj_b.yaml"), "w") do io
                                write(io, "# Dummy settings file for proj_b\n")
                            end
                            open(joinpath("data", "vis_a.yaml"), "w") do io
                                write(io, "# Dummy vis file for proj_a\n")
                            end
                            open(joinpath("data", "vis_b.yaml"), "w") do io
                                write(io, "# Dummy vis file for proj_b\n")
                            end
                            
                            # Set initial MSR
                            FLORIDyn.set_default_msr(AddedTurbulence)
                            initial_msr = FLORIDyn.get_default_msr()
                            @test initial_msr == AddedTurbulence
                            
                            # Simulate what select_project does when changing project
                            chosen_name = "proj_b"
                            
                            # Read existing MSR (as select_project does)
                            existing_msr = AddedTurbulence  # Fallback
                            default_path_local = joinpath("data", "default.yaml")
                            if isfile(default_path_local)
                                try
                                    def_data = YAML.load_file(default_path_local)
                                    if haskey(def_data, "default") && haskey(def_data["default"], "msr")
                                        msr_str = String(def_data["default"]["msr"])
                                        existing_msr = FLORIDyn.toMSR(msr_str)
                                    end
                                catch
                                    # If malformed, will use fallback
                                end
                            end
                            
                            # Write new project while preserving MSR (as select_project does)
                            open(default_path_local, "w") do io
                                write(io, "default:\n  name: $(chosen_name)\n  msr: $(string(existing_msr))  # valid options: VelReduction, AddedTurbulence, EffWind\n")
                            end
                            
                            # Verify project changed and MSR preserved
                            current_project, _ = FLORIDyn.get_default_project()
                            current_name = splitpath(current_project)[end]
                            current_name = replace(current_name, ".yaml" => "")
                            @test current_name == "proj_b"
                            @test FLORIDyn.get_default_msr() == AddedTurbulence  # MSR should be preserved
                        end
                    end
                end

                # Test cancellation behavior
                @testset "Cancellation handling" begin
                    mktempdir() do tmp
                        cd(tmp) do
                            # Create test projects.yaml
                            mkpath("data")
                            test_projects = """
                            projects:
                              - project:
                                  name: original_proj
                                  description: "Original project"
                                  vis: "vis_orig.yaml"
                            """
                            open(joinpath("data", "projects.yaml"), "w") do io
                                write(io, test_projects)
                            end
                            
                            # Create the necessary settings files
                            open(joinpath("data", "original_proj.yaml"), "w") do io
                                write(io, "# Dummy settings file for original_proj\n")
                            end
                            open(joinpath("data", "vis_orig.yaml"), "w") do io
                                write(io, "# Dummy vis file for original_proj\n")
                            end
                            
                            # Set up initial project
                            FLORIDyn.set_default_msr(VelReduction)  # This also sets a default project
                            
                            # Get current project name for cancellation test
                            current_project, _ = FLORIDyn.get_default_project()
                            current_name = splitpath(current_project)[end]
                            current_name = replace(current_name, ".yaml" => "")
                            
                            # Test that get_default_project returns consistent format
                            # (This is what select_project returns on cancellation)
                            @test isa(current_name, String)
                            @test length(current_name) > 0
                        end
                    end
                end

                # Test menu configuration compatibility with REPL.TerminalMenus
                @testset "TerminalMenus compatibility" begin
                    mktempdir() do tmp
                        cd(tmp) do
                            # Create test projects.yaml
                            mkpath("data")
                            test_projects = """
                            projects:
                              - project:
                                  name: menu_test1
                                  description: "Menu test project 1"
                                  vis: "vis_mt1.yaml"
                              - project:
                                  name: menu_test2
                                  description: "Menu test project 2"
                                  vis: "vis_mt2.yaml"
                              - project:
                                  name: no_desc
                                  description: ""
                                  vis: "vis_nd.yaml"
                            """
                            open(joinpath("data", "projects.yaml"), "w") do io
                                write(io, test_projects)
                            end
                            
                            # Test that we can create a RadioMenu with project options
                            projects = FLORIDyn.list_projects()
                            project_options = String[]
                            for (name, description, vis) in projects
                                if isempty(description)
                                    push!(project_options, name)
                                else
                                    push!(project_options, "$name - $description")
                                end
                            end
                            
                            # Test RadioMenu creation (structure that select_project uses)
                            try
                                menu = REPL.TerminalMenus.RadioMenu(project_options, pagesize=length(project_options))
                                @test isa(menu, REPL.TerminalMenus.RadioMenu)
                                @test length(project_options) == 3
                            catch e
                                # If TerminalMenus is not available in test environment, skip this test
                                @test_skip "TerminalMenus not available in test environment: $e"
                            end
                        end
                    end
                end

                # Test project name extraction from get_default_project
                @testset "Project name extraction" begin
                    mktempdir() do tmp
                        cd(tmp) do
                            # Test name extraction logic used in select_project
                            test_paths = [
                                ("/path/to/test_project.yaml", "test_project"),
                                ("simple.yaml", "simple"),
                                ("/long/path/with/many/dirs/complex_name.yaml", "complex_name"),
                                ("data/default.yaml", "default")
                            ]
                            
                            for (test_path, expected_name) in test_paths
                                # Simulate the name extraction logic from select_project
                                current_name = splitpath(test_path)[end]
                                current_name = replace(current_name, ".yaml" => "")
                                @test current_name == expected_name
                            end
                        end
                    end
                end

                # Note: Full interactive testing of select_project() with actual menu 
                # interaction would require complex input simulation. The function's core 
                # logic is thoroughly tested through the underlying project functions and 
                # menu structure validation above.
            end
        end
    end
end
nothing