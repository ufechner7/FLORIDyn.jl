# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Test
using FLORIDyn
using Pkg

@testset "Copy Functions Tests" begin
    # Create a temporary directory for testing
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        @testset "copy_bin function" begin
            # Test that copy_bin creates bin directory and copies run_julia script
            @test !isdir("bin")
            
            # Call copy_bin function (not exported, so call directly)
            FLORIDyn.copy_bin()
            
            # Check that bin directory was created
            @test isdir("bin")
            
            # Check that run_julia script was copied
            @test isfile("bin/run_julia")
            
            # Check that the file has executable permissions
            stat_info = stat("bin/run_julia")
            @test (stat_info.mode & 0o111) != 0  # Check if any execute bits are set
            
            # Check that the file content is copied correctly
            # Read the original file from the package
            pkg_path = pkgdir(FLORIDyn)
            original_file = joinpath(pkg_path, "bin", "run_julia")
            if isfile(original_file)
                original_content = read(original_file, String)
                copied_content = read("bin/run_julia", String)
                @test original_content == copied_content
            end
            
            # Test that calling copy_bin again doesn't fail (should overwrite)
            @test_nowarn FLORIDyn.copy_bin()
            @test isfile("bin/run_julia")
        end
        
        @testset "copy_bin error handling" begin
            # Test behavior when source package path cannot be found
            # This is hard to test directly, but we can at least ensure
            # the function doesn't crash unexpectedly
            @test_nowarn FLORIDyn.copy_bin()
        end
        
        @testset "copy_bin with existing bin directory" begin
            # Test when bin directory already exists
            mkpath("bin2")
            cd("bin2")
            mkdir("bin")
            
            # Should still work and copy the file
            @test_nowarn FLORIDyn.copy_bin()
            @test isfile("bin/run_julia")
        end
        
    finally
        cd(original_dir)
        # Clean up: remove the temporary directory
        rm(test_dir, recursive=true, force=true)
    end
end

@testset "copy_bin permissions test" begin
    # Test file permissions more thoroughly
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        FLORIDyn.copy_bin()
        
        if isfile("bin/run_julia")
            stat_info = stat("bin/run_julia")
            # Check that the file has the expected permissions (0o774)
            # On some systems, umask might affect this, so we check for reasonable permissions
            mode = stat_info.mode & 0o777
            
            # Should have read and execute for owner and group, at minimum
            @test (mode & 0o500) == 0o500  # Owner read and execute
            @test (mode & 0o050) == 0o050  # Group read and execute
            
            # Should be executable by owner
            @test (mode & 0o100) != 0
        end
        
    finally
        cd(original_dir)
        rm(test_dir, recursive=true, force=true)
    end
end

@testset "copy_bin integration test" begin
    # Test that the copied script is actually functional
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        FLORIDyn.copy_bin()
        
        if isfile("bin/run_julia")
            # Test that the script exists and has a shebang or julia reference
            content = read("bin/run_julia", String)
            @test occursin("julia", lowercase(content)) || occursin("#!/", content)
        end
        
    finally
        cd(original_dir)
        rm(test_dir, recursive=true, force=true)
    end
end

@testset "copy_model_settings function" begin
    # Create a temporary directory for testing
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        @testset "basic functionality" begin
            # Test that copy_model_settings creates data directory and copies files
            @test !isdir("data")
            
            # Call copy_model_settings function (not exported, so call directly)
            FLORIDyn.copy_model_settings()
            
            # Check that data directory was created
            @test isdir("data")
            
            # Check that the YAML configuration file was copied
            @test isfile("data/2021_9T_Data.yaml")
            
            # Check that the 2021_9T_Data directory was copied
            @test isdir("data/2021_9T_Data")
            
            # Check that the file has correct permissions
            stat_info = stat("data/2021_9T_Data.yaml")
            mode = stat_info.mode & 0o777
            @test (mode & 0o700) == 0o700  # Owner should have rwx
        end
        
        @testset "content validation" begin
            # Call the function first
            FLORIDyn.copy_model_settings()
            
            # Check that the YAML file content is copied correctly
            pkg_path = pkgdir(FLORIDyn)
            original_yaml = joinpath(pkg_path, "data", "2021_9T_Data.yaml")
            if isfile(original_yaml)
                original_content = read(original_yaml, String)
                copied_content = read("data/2021_9T_Data.yaml", String)
                @test original_content == copied_content
            end
            
            # Check that 2021_9T_Data directory contains expected files
            data_dir = "data/2021_9T_Data"
            @test isdir(data_dir)
            
            # Check for some expected CSV files (these should exist based on the function documentation)
            expected_files = ["SOWFA_bladePitch.csv", "SOWFA_generatorPower.csv", "U.csv", "WindDir.csv"]
            for file in expected_files
                file_path = joinpath(data_dir, file)
                if isfile(joinpath(pkg_path, "data", "2021_9T_Data", file))
                    @test isfile(file_path)
                end
            end
        end
        
        @testset "overwrite behavior" begin
            # First call
            FLORIDyn.copy_model_settings()
            @test isfile("data/2021_9T_Data.yaml")
            
            # Modify the copied file
            open("data/2021_9T_Data.yaml", "a") do f
                write(f, "\n# Modified for testing")
            end
            original_size = filesize("data/2021_9T_Data.yaml")
            
            # Second call should overwrite
            FLORIDyn.copy_model_settings()
            @test isfile("data/2021_9T_Data.yaml")
            
            # File should be restored to original content (smaller size)
            new_size = filesize("data/2021_9T_Data.yaml")
            @test new_size < original_size  # Should be back to original
        end
        
        @testset "directory permissions" begin
            FLORIDyn.copy_model_settings()
            
            # Check permissions on files in the 2021_9T_Data directory
            data_dir = "data/2021_9T_Data"
            if isdir(data_dir)
                for (root, dirs, files_in_dir) in walkdir(data_dir)
                    for file in files_in_dir
                        file_path = joinpath(root, file)
                        stat_info = stat(file_path)
                        mode = stat_info.mode & 0o777
                        # Files should have 0o774 permissions as set by the function
                        @test (mode & 0o700) == 0o700  # Owner should have rwx
                        @test (mode & 0o070) == 0o070  # Group should have rwx
                        @test (mode & 0o004) == 0o004  # Others should have read
                    end
                end
            end
        end
        
        @testset "error handling with existing data directory" begin
            # Create data directory first
            mkpath("data")
            existing_file = "data/existing_file.txt"
            write(existing_file, "This file existed before")
            
            # Should still work and copy the model settings
            @test_nowarn FLORIDyn.copy_model_settings()
            @test isfile("data/2021_9T_Data.yaml")
            @test isdir("data/2021_9T_Data")
            @test isfile(existing_file)  # Existing file should remain
        end
        
        @testset "function output verification" begin
            # Test that the function runs without error and produces output
            # Note: The actual output verification is hard to test reliably across systems
            @test_nowarn FLORIDyn.copy_model_settings()
            @test isfile("data/2021_9T_Data.yaml")
            @test isdir("data/2021_9T_Data")
        end
        
    finally
        cd(original_dir)
        # Clean up: remove the temporary directory
        rm(test_dir, recursive=true, force=true)
    end
end

@testset "copy_model_settings integration test" begin
    # Test that copied files are valid and can be used
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        FLORIDyn.copy_model_settings()
        
        # Test that the YAML file is valid
        if isfile("data/2021_9T_Data.yaml")
            # Try to read the YAML file (basic validation)
            content = read("data/2021_9T_Data.yaml", String)
            @test !isempty(content)
            @test occursin("yaml", lowercase(content)) || occursin("wind", lowercase(content)) || length(content) > 100
        end
        
        # Test that the copied data directory structure is reasonable
        data_dir = "data/2021_9T_Data"
        if isdir(data_dir)
            files_in_dir = readdir(data_dir)
            @test length(files_in_dir) > 0  # Should contain some files
            
            # Check that at least some CSV files exist
            csv_files = filter(f -> endswith(f, ".csv"), files_in_dir)
            @test length(csv_files) > 0  # Should have some CSV files
        end
        
    finally
        cd(original_dir)
        rm(test_dir, recursive=true, force=true)
    end
end

@testset "copy_examples function" begin
    # Create a temporary directory for testing
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        @testset "basic functionality" begin
            # Test that copy_examples creates examples directory and copies files
            @test !isdir("examples")
            
            # Call copy_examples function (not exported, so call directly)
            FLORIDyn.copy_examples()
            
            # Check that examples directory was created
            @test isdir("examples")
            
            # Check that the expected example files were copied
            expected_files = ["main.jl", "menu.jl", "create_videos.jl", "play_video.jl"]
            for file in expected_files
                @test isfile(joinpath("examples", file))
            end
            
            # Check that files have correct permissions (0o774)
            for file in expected_files
                if isfile(joinpath("examples", file))
                    stat_info = stat(joinpath("examples", file))
                    mode = stat_info.mode & 0o777
                    @test (mode & 0o700) == 0o700  # Owner should have rwx
                    @test (mode & 0o070) == 0o070  # Group should have rwx
                    @test (mode & 0o004) == 0o004  # Others should have read
                end
            end
        end
        
        @testset "content validation" begin
            # Call the function first
            FLORIDyn.copy_examples()
            
            # Check that the Julia files content is copied correctly
            pkg_path = pkgdir(FLORIDyn)
            examples_files = ["main.jl", "menu.jl", "create_videos.jl", "play_video.jl"]
            
            for file in examples_files
                original_file = joinpath(pkg_path, "examples", file)
                copied_file = joinpath("examples", file)
                
                if isfile(original_file) && isfile(copied_file)
                    original_content = read(original_file, String)
                    copied_content = read(copied_file, String)
                    @test original_content == copied_content
                end
            end
            
            # Check that the files are valid Julia files
            for file in examples_files
                file_path = joinpath("examples", file)
                if isfile(file_path)
                    content = read(file_path, String)
                    @test !isempty(content)
                    # Basic check for Julia syntax
                    @test occursin("using", content) || occursin("function", content) || occursin("# ", content)
                end
            end
        end
        
        @testset "overwrite behavior" begin
            # First call
            FLORIDyn.copy_examples()
            @test isfile("examples/main.jl")
            
            # Modify one of the copied files
            main_file = "examples/main.jl"
            open(main_file, "a") do f
                write(f, "\n# Modified for testing")
            end
            original_size = filesize(main_file)
            
            # Second call should overwrite
            FLORIDyn.copy_examples()
            @test isfile(main_file)
            
            # File should be restored to original content (smaller size)
            new_size = filesize(main_file)
            @test new_size < original_size
        end
        
        @testset "directory permissions and structure" begin
            FLORIDyn.copy_examples()
            
            # Check that the examples directory exists and is accessible
            @test isdir("examples")
            @test isreadable("examples")
            
            # Check that we can list files in the directory
            files_in_dir = readdir("examples")
            @test length(files_in_dir) > 0
            
            # Check that all files are readable
            for file in files_in_dir
                file_path = joinpath("examples", file)
                if isfile(file_path)
                    @test isreadable(file_path)
                    
                    # Check permissions
                    stat_info = stat(file_path)
                    mode = stat_info.mode & 0o777
                    @test (mode & 0o004) != 0
                end
            end
        end
        
        @testset "error handling with existing examples directory" begin
            # Create examples directory first with an existing file
            mkpath("examples")
            existing_file = "examples/existing_file.jl"
            write(existing_file, "# This file existed before")
            
            # Should still work and copy the example files
            @test_nowarn FLORIDyn.copy_examples()
            @test isfile("examples/main.jl")
            @test isfile("examples/menu.jl")
            @test isfile(existing_file)
        end
        
        @testset "function robustness" begin
            # Test multiple calls don't cause issues
            @test_nowarn FLORIDyn.copy_examples()
            @test_nowarn FLORIDyn.copy_examples()
            @test_nowarn FLORIDyn.copy_examples()
            
            # Check that files are still there and valid
            @test isdir("examples")
            @test isfile("examples/main.jl")
            @test isfile("examples/menu.jl")
        end
        
        @testset "file count consistency" begin
            FLORIDyn.copy_examples()
            
            # Get list of files from package examples directory
            pkg_path = pkgdir(FLORIDyn)
            src_examples = joinpath(pkg_path, "examples")
            
            if isdir(src_examples)
                src_files = readdir(src_examples)
                copied_files = readdir("examples")
                
                # Filter out any test files we might have created during testing
                julia_src_files = filter(f -> endswith(f, ".jl"), src_files)
                julia_copied_files = filter(f -> endswith(f, ".jl") && !startswith(f, "existing_"), copied_files)
                
                # Should have copied all Julia files from source
                @test length(julia_copied_files) >= length(julia_src_files)
                
                # Check that all source files are present in copied directory
                for src_file in julia_src_files
                    @test src_file in copied_files
                end
            end
        end
        
    finally
        cd(original_dir)
        # Clean up: remove the temporary directory
        rm(test_dir, recursive=true, force=true)
    end
end

@testset "copy_examples integration test" begin
    # Test that copied example files are valid and can be used
    test_dir = mktempdir()
    original_dir = pwd()
    
    try
        cd(test_dir)
        
        FLORIDyn.copy_examples()
        
        # Test that the Julia files are syntactically valid
        examples_files = ["main.jl", "menu.jl", "create_videos.jl", "play_video.jl"]
        
        for file in examples_files
            file_path = joinpath("examples", file)
            if isfile(file_path)
                content = read(file_path, String)
                @test !isempty(content)
                
                # Check for common Julia patterns
                julia_patterns = ["using", "function", "if", "end", "#"]
                has_julia_pattern = any(pattern -> occursin(pattern, content), julia_patterns)
                @test has_julia_pattern
                
                # Try to parse the file as Julia code (basic syntax check)
                # This is a more thorough test but might be sensitive to syntax
                try
                    # We can't easily parse without potentially executing, so we do basic checks
                    @test !occursin("<<<<<<< HEAD", content)
                    @test !occursin(">>>>>>> ", content)
                catch e
                    @warn "Could not perform syntax check on $file: $e"
                end
            end
        end
        
        # Test that main.jl contains expected content
        main_file = "examples/main.jl"
        if isfile(main_file)
            content = read(main_file, String)
            @test occursin("FLORIDyn", content)
            # Check for typical main.jl patterns
            @test occursin("using", content) || occursin("import", content)
        end
        
        # Test that the examples directory structure is reasonable
        @test isdir("examples")
        files_in_dir = readdir("examples")
        @test length(files_in_dir) >= 3
        
        # Check that files have reasonable sizes (not empty, not too large)
        for file in files_in_dir
            file_path = joinpath("examples", file)
            if isfile(file_path) && endswith(file, ".jl")
                size = filesize(file_path)
                @test size > 50
                @test size < 100000
            end
        end
        
    finally
        cd(original_dir)
        rm(test_dir, recursive=true, force=true)
    end
end
