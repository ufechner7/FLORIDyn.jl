# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test
using Dates

@testset verbose=true "High-Resolution Time Functions" begin
    
    @testset "now_microseconds" begin
        # Test that the function returns a string
        result = now_microseconds()
        @test isa(result, String)
        
        # Test the format pattern: YYYY-mm-ddTHH-MM-SS.uuuuuu
        @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{6}$", result)
        
        # Test that successive calls return different microseconds (most of the time)
        result1 = now_microseconds()
        sleep(0.001)  # Sleep 1ms to ensure different microseconds
        result2 = now_microseconds()
        @test result1 != result2
        
        # Test that the timestamp is reasonable (within a few seconds of now)
        timestamp_part = split(result, ".")[1]  # Get the part before microseconds
        # Convert hyphens to colons in time portion for DateTime parsing
        datetime_str = replace(timestamp_part, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
        parsed_time = DateTime(datetime_str, "yyyy-mm-ddTHH:MM:SS")
        current_time = now()
        # Calculate time difference in seconds using Millisecond conversion
        time_diff_ms = abs(Millisecond(current_time - parsed_time).value)
        time_diff = time_diff_ms / 1000.0  # Convert to seconds
        @test time_diff < 5.0  # Should be within 5 seconds
        
        # Test that microseconds part has exactly 6 digits
        microseconds_part = split(result, ".")[2]
        @test length(microseconds_part) == 6
        @test all(isdigit, microseconds_part)
    end

    @testset "now_nanoseconds" begin
        # Test that the function returns a string
        result = now_nanoseconds()
        @test isa(result, String)
        
        # Test the format pattern: YYYY-mm-ddTHH-MM-SS.nnnnnnnnn
        @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{9}$", result)
        
        # Test that successive calls return different nanoseconds (most of the time)
        result1 = now_nanoseconds()
        result2 = now_nanoseconds()
        @test result1 != result2  # Should be different due to nanosecond precision
        
        # Test that nanoseconds part has exactly 9 digits
        nanoseconds_part = split(result, ".")[2]
        @test length(nanoseconds_part) == 9
        @test all(isdigit, nanoseconds_part)
        
        # Test that the timestamp is reasonable (within a few seconds of now)
        timestamp_part = split(result, ".")[1]  # Get the part before nanoseconds
        datetime_str = replace(timestamp_part, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
        parsed_time = DateTime(datetime_str, "yyyy-mm-ddTHH:MM:SS")
        current_time = now()
        time_diff_ms = abs(Millisecond(current_time - parsed_time).value)
        time_diff = time_diff_ms / 1000.0
        @test time_diff < 5.0
        
        # Test that nanosecond precision is actually higher than microsecond
        micro_result = now_microseconds()
        nano_result = now_nanoseconds()
        # The nanosecond version should have 3 more digits of precision
        @test length(split(nano_result, ".")[2]) == length(split(micro_result, ".")[2]) + 3
    end
    
    @testset "precise_now" begin
        # Test that the function returns a String 
        result = precise_now()
        @test isa(result, String)
        
        # Test that the format matches expected pattern (YYYY-MM-DDTHH-MM-SS.microseconds)
        @test match(r"^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{6}$", result) !== nothing
        
        # Test consecutive calls to verify precision
        result1 = precise_now()
        result2 = precise_now()
        @test result1 <= result2  # Should be monotonic or equal (lexicographically)
    end

    @testset "unique_name" begin
        # Test basic functionality
        name1 = unique_name()
        @test isa(name1, String)
        @test startswith(name1, "floridyn_run_")
        
        # Test uniqueness
        name2 = unique_name()
        @test name1 != name2
        
        # Test format - should contain timestamp
        @test occursin(r"^floridyn_run_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{6}$", name1)
        @test occursin(r"^floridyn_run_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{6}$", name2)
        
        # Test that multiple rapid calls generate unique names
        names = [unique_name() for _ in 1:5]
        @test length(unique(names)) == 5  # All names should be unique
        
        # Test that all names have the correct prefix
        for name in names
            @test startswith(name, "floridyn_run_")
        end
        
        # Test timestamp parsing - names should be parseable timestamps
        for name in names
            timestamp_str = replace(name, "floridyn_run_" => "")
            timestamp_part = split(timestamp_str, ".")[1]
            datetime_str = replace(timestamp_part, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
            @test_nowarn DateTime(datetime_str, "yyyy-mm-ddTHH:MM:SS")
        end
    end
    
    @testset "Time precision and uniqueness" begin
        # Test microsecond precision guarantees uniqueness
        timestamps = [now_microseconds() for _ in 1:10]
        @test length(unique(timestamps)) == length(timestamps)
        
        # Test nanosecond precision guarantees uniqueness  
        nano_timestamps = [now_nanoseconds() for _ in 1:10]
        @test length(unique(nano_timestamps)) == length(nano_timestamps)
        
        # Test unique_name uniqueness even with rapid calls
        folder_names = [unique_name() for _ in 1:10]
        @test length(unique(folder_names)) == length(folder_names)
    end
    
    @testset "Format consistency" begin
        # Test that all timestamp functions use consistent date formatting
        micro_time = now_microseconds()
        nano_time = now_nanoseconds()
        unique_folder = unique_name()
        
        # Extract date parts and verify they're the same (within same second)
        micro_date = split(micro_time, "T")[1]
        nano_date = split(nano_time, "T")[1] 
        unique_date = split(replace(unique_folder, "floridyn_run_" => ""), "T")[1]
        
        @test micro_date == nano_date
        # Unique date should be same or very close (might differ by seconds if crossing boundary)
        @test abs(DateTime(unique_date) - DateTime(micro_date)) < Millisecond(2000)
        
        # Test time format consistency (HH-MM-SS pattern)
        micro_time_part = split(split(micro_time, "T")[2], ".")[1]
        nano_time_part = split(split(nano_time, "T")[2], ".")[1]
        @test occursin(r"^\d{2}-\d{2}-\d{2}$", micro_time_part)
        @test occursin(r"^\d{2}-\d{2}-\d{2}$", nano_time_part)
    end

    @testset "Edge cases and robustness" begin
        # Test multiple rapid successive calls
        results = []
        for _ in 1:100
            push!(results, now_microseconds())
        end
        @test length(unique(results)) == length(results)  # All should be unique
        
        # Test that functions work across second boundaries
        # (This is probabilistic, but should work most of the time)
        start_time = now()
        results = []
        while Millisecond(now() - start_time).value < 2000  # Run for 2 seconds
            push!(results, now_microseconds())
            if length(results) > 50  # Don't run forever
                break
            end
        end
        @test length(unique(results)) == length(results)
        
        # Test format robustness
        for result in results[1:min(10, length(results))]
            @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{6}$", result)
        end
    end

    @testset "delete_results function" begin
        @testset "Input validation" begin
            # Create a test Vis object with unique identifiers
            unique_id = now_nanoseconds()
            vis = Vis(online=false, output_folder="test_validation_$(unique_id)", video_folder="test_validation_video_$(unique_id)")
            
            # Test with non-positive n
            @test_logs (:warn, r"Number of folders to delete must be positive") begin
                result = delete_results(vis, 0)
                @test isempty(result)
            end
            
            @test_logs (:warn, r"Number of folders to delete must be positive") begin
                result = delete_results(vis, -1)
                @test isempty(result)
            end
        end
        
        @testset "Empty directory behavior" begin
            # Create a test Vis object with isolated directories
            unique_id = rand(UInt32)
            vis = Vis(online=false, output_folder="test_empty_output_$(unique_id)", video_folder="test_empty_video_$(unique_id)")
            
            # Test with empty directory (no floridyn_run folders)
            @test_logs (:info, r"Directory does not exist") match_mode=:any begin
                result = delete_results(vis, 1)
                @test isempty(result)
            end
            
            # Test dry run with empty directory
            @test_logs (:info, r"Directory does not exist") match_mode=:any begin
                result = delete_results(vis, 1, dry_run=true)
                @test isempty(result)
            end
        end

        @testset "Single folder operations" begin
            # Create a test Vis object with isolated paths (use unique names for test isolation)
            unique_id = rand(UInt32)
            vis = Vis(online=false, output_folder="test_single_$(unique_id)", video_folder="test_single_video_$(unique_id)")  
            
            # Get both paths - function now uses base directories directly
            output_dir = joinpath(pwd(), vis.output_folder)
            video_dir = joinpath(pwd(), vis.video_folder)
            
            # Create the directory structures manually
            mkpath(output_dir)
            test_folder = joinpath(output_dir, "floridyn_run_2025-01-01T12-00-00.123456")
            mkdir(test_folder)
            
            # Create a dummy file inside to make it non-empty
            touch(joinpath(test_folder, "test_file.txt"))
            
            # Test dry run - should find the folder but not delete it
            @test_logs (:info, r"DRY RUN") match_mode=:any begin
                result = delete_results(vis, 1, dry_run=true)
                @test length(result) == 1
                @test basename(result[1]) == "floridyn_run_2025-01-01T12-00-00.123456"
                @test isdir(test_folder)  # Should still exist
            end
            
            # Test actual deletion
            @test_logs (:info, r"Deleting") match_mode=:any begin
                result = delete_results(vis, 1)
                @test length(result) == 1
                @test basename(result[1]) == "floridyn_run_2025-01-01T12-00-00.123456"
                @test !isdir(test_folder)  # Should be deleted
            end
            
            # Clean up both directories
            rm(joinpath(pwd(), vis.output_folder), recursive=true, force=true)
            rm(joinpath(pwd(), vis.video_folder), recursive=true, force=true)
        end
            
        @testset "Multiple folder operations" begin
            # Create a test Vis object with isolated paths
            unique_id = rand(UInt32)
            vis = Vis(online=false, output_folder="test_multiple_$(unique_id)", video_folder="test_multiple_video_$(unique_id)")
            
            # Get both paths - function now uses base directories directly
            output_dir = joinpath(pwd(), vis.output_folder)
            video_dir = joinpath(pwd(), vis.video_folder)
            mkpath(output_dir)
            
            # Create multiple test directories with different timestamps in output_dir
            test_folders = [
                "floridyn_run_2025-01-01T10-00-00.111111",  # Oldest
                "floridyn_run_2025-01-01T11-00-00.222222",  # Middle
                "floridyn_run_2025-01-01T12-00-00.333333",  # Newest
            ]
            
            folder_paths = String[]
            for folder_name in test_folders
                folder_path = joinpath(output_dir, folder_name)
                mkdir(folder_path)
                touch(joinpath(folder_path, "test_file.txt"))
                push!(folder_paths, folder_path)
                # Small delay to ensure different modification times
                sleep(0.01)
            end
            
            # Test dry run for deleting 2 newest
            @test_logs (:info, r"DRY RUN") match_mode=:any begin
                result = delete_results(vis, 2, dry_run=true)
                @test length(result) >= 2  # Should find at least 2
                # All folders should still exist
                @test all(isdir, folder_paths)
            end
            
            # Test actual deletion of 2 newest
            @test_logs (:info, r"Deleting") match_mode=:any begin
                result = delete_results(vis, 2)
                @test length(result) >= 2  # Should delete at least 2
            end
            
            # Clean up both directories
            rm(joinpath(pwd(), vis.output_folder), recursive=true, force=true)
            rm(joinpath(pwd(), vis.video_folder), recursive=true, force=true)
        end
        
        @testset "Boundary conditions" begin
            # Create a test Vis object
            unique_id = rand(UInt32)
            vis = Vis(online=false, output_folder="test_boundary_$(unique_id)", video_folder="test_boundary_video_$(unique_id)")
            
            # Get the paths and create directory structure
            output_dir = joinpath(pwd(), vis.output_folder)
            mkpath(output_dir)
            
            # Create 3 test directories
            test_folders = [
                "floridyn_run_2025-01-01T10-00-00.111111",
                "floridyn_run_2025-01-01T11-00-00.222222", 
                "floridyn_run_2025-01-01T12-00-00.333333",
            ]
            
            folder_paths = String[]
            for folder_name in test_folders
                folder_path = joinpath(output_dir, folder_name)
                mkdir(folder_path)
                push!(folder_paths, folder_path)
                sleep(0.01)
            end
            
            # Test requesting more folders than available
            @test_logs (:info, r"Deleting") match_mode=:any begin
                result = delete_results(vis, 10)  # Request 10, only 3 available in output
                @test length(result) >= 3  # Should delete at least the 3 from output
            end
            
            # Clean up
            rm(joinpath(pwd(), vis.output_folder), recursive=true, force=true)
        end
        
        @testset "Error handling" begin
            # Create a test Vis object
            unique_id = rand(UInt32)
            vis = Vis(online=false, output_folder="test_error_$(unique_id)", video_folder="test_error_video_$(unique_id)")
            
            # Test deletion - should handle non-existent directories gracefully and return empty results
            result = delete_results(vis, 1)
            @test result == String[]  # Should return empty array when no directories exist
            @test isa(result, Vector{String})  # Should return the correct type
        end
        
        @testset "Integration with find_floridyn_runs" begin
            # Create a test Vis object
            unique_id = rand(UInt32)
            vis = Vis(online=false, output_folder="test_integration_$(unique_id)", video_folder="test_integration_video_$(unique_id)")
            
            # Get the paths and create directory structure
            output_dir = joinpath(pwd(), vis.output_folder)
            mkpath(output_dir)
            
            # Create mixed directories (some floridyn_run, some not)
            mixed_dirs = [
                "floridyn_run_2025-01-01T10-00-00.111111",
                "other_folder", 
                "floridyn_run_2025-01-01T11-00-00.222222",
                "not_floridyn_folder",
                "floridyn_run_2025-01-01T12-00-00.333333"
            ]
            
            for dir_name in mixed_dirs
                mkdir(joinpath(output_dir, dir_name))
                sleep(0.01)
            end
            
            # Should only find and delete floridyn_run directories
            result = delete_results(vis, 10)
            @test length(result) >= 3  # Should find at least the 3 floridyn_run directories
            
            # Check that non-floridyn directories are preserved
            @test isdir(joinpath(output_dir, "other_folder"))
            @test isdir(joinpath(output_dir, "not_floridyn_folder"))
            
            # Check that floridyn_run directories are deleted
            @test !isdir(joinpath(output_dir, "floridyn_run_2025-01-01T10-00-00.111111"))
            @test !isdir(joinpath(output_dir, "floridyn_run_2025-01-01T11-00-00.222222"))
            @test !isdir(joinpath(output_dir, "floridyn_run_2025-01-01T12-00-00.333333"))
            
            # Clean up
            rm(joinpath(pwd(), vis.output_folder), recursive=true, force=true)
        end
        
        @testset "Default parameter behavior" begin
            # Create a test Vis object
            unique_id = rand(UInt32)
            vis = Vis(online=false, output_folder="test_default_$(unique_id)", video_folder="test_default_video_$(unique_id)")
            
            # Get the paths and create directory structure
            output_dir = joinpath(pwd(), vis.output_folder)
            mkpath(output_dir)
            
            # Test default n=1
            test_folder = joinpath(output_dir, "floridyn_run_2025-01-01T12-00-00.123456")
            mkdir(test_folder)
            
            # Test that calling without explicit n deletes 1 folder
            result = delete_results(vis)  # Uses default n=1
            @test length(result) >= 1
            @test !isdir(test_folder)
            
            # Clean up
            rm(joinpath(pwd(), vis.output_folder), recursive=true, force=true)
        end
    end
end
