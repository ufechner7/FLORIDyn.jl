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
        # Parse the timestamp to check it's close to current time
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
        # No sleep needed for nanosecond precision
        result2 = now_nanoseconds()
        @test result1 != result2  # Should be different due to nanosecond precision
        
        # Test that the timestamp is reasonable (within a few seconds of now)
        timestamp_part = split(result, ".")[1]  # Get the part before nanoseconds
        # Convert hyphens to colons in time portion for DateTime parsing
        datetime_str = replace(timestamp_part, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
        parsed_time = DateTime(datetime_str, "yyyy-mm-ddTHH:MM:SS")
        current_time = now()
        # Calculate time difference in seconds using Millisecond conversion
        time_diff_ms = abs(Millisecond(current_time - parsed_time).value)
        time_diff = time_diff_ms / 1000.0  # Convert to seconds
        @test time_diff < 5.0  # Should be within 5 seconds
        
        # Test that nanoseconds part has exactly 9 digits
        nanoseconds_part = split(result, ".")[2]
        @test length(nanoseconds_part) == 9
        @test all(isdigit, nanoseconds_part)
    end
    
    @testset "precise_now" begin
        # Test that precise_now is equivalent to now_microseconds
        result_precise = precise_now()
        result_micro = now_microseconds()
        
        @test isa(result_precise, String)
        
        # Test same format as now_microseconds
        @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{6}$", result_precise)
        
        # Test that both functions produce similar timestamps (within a small time window)
        # Get the microsecond parts to compare
        precise_microseconds = parse(Int, split(result_precise, ".")[2])
        micro_microseconds = parse(Int, split(result_micro, ".")[2])
        
        # The difference should be small (both should be very close in time)
        # Allow for some variance due to execution time
        @test abs(precise_microseconds - micro_microseconds) < 100000  # Within 0.1 second
    end
    
    @testset "unique_name" begin
        # Test that the function returns a string with correct prefix
        result = unique_name()
        @test isa(result, String)
        @test startswith(result, "floridyn_run_")
        
        # Test the format pattern: floridyn_run_YYYY-mm-ddTHH-MM-SS.uuuuuu
        @test occursin(r"^floridyn_run_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{6}$", result)
        
        # Test that successive calls return unique names
        result1 = unique_name()
        sleep(0.001)  # Small delay to ensure different timestamps
        result2 = unique_name()
        @test result1 != result2
        
        # Test that multiple calls in quick succession produce unique names
        names = Set{String}()
        for i in 1:10
            push!(names, unique_name())
        end
        @test length(names) == 10  # All names should be unique
        
        # Test that the timestamp part is reasonable
        timestamp_part = result[14:end]  # Remove "floridyn_run_" prefix (13 chars + 1)
        timestamp_base = split(timestamp_part, ".")[1]
        # Convert hyphens to colons in time portion for DateTime parsing
        datetime_str = replace(timestamp_base, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
        parsed_time = DateTime(datetime_str, "yyyy-mm-ddTHH:MM:SS")
        current_time = now()
        # Calculate time difference in seconds using Millisecond conversion
        time_diff_ms = abs(Millisecond(current_time - parsed_time).value)
        time_diff = time_diff_ms / 1000.0  # Convert to seconds
        @test time_diff < 5.0  # Should be within 5 seconds
        
        # Test that names are suitable as directory names (no invalid characters)
        @test !occursin(r"[<>:\"/\\|?*]", result)  # No Windows-invalid chars
        @test !occursin(" ", result)  # No spaces
    end
    
    @testset "Time precision and uniqueness" begin
        # Test that microsecond precision is sufficient for rapid calls
        microsecond_times = [now_microseconds() for _ in 1:100]
        @test length(unique(microsecond_times)) >= 90  # Most should be unique
        
        # Test that nanosecond precision gives better uniqueness
        nanosecond_times = [now_nanoseconds() for _ in 1:100]
        @test length(unique(nanosecond_times)) >= 95  # Most should be unique
        
        # Test that unique_name generates unique names even in rapid succession
        unique_names = [unique_name() for _ in 1:50]
        @test length(unique(unique_names)) == 50  # All should be unique
    end
    
    @testset "Format consistency" begin
        # Test that all functions produce consistent date formatting
        micro_result = now_microseconds()
        nano_result = now_nanoseconds()
        precise_result = precise_now()
        unique_result = unique_name()
        
        # Extract the base timestamp parts (before fractional seconds)
        micro_base = split(micro_result, ".")[1]
        nano_base = split(nano_result, ".")[1]
        precise_base = split(precise_result, ".")[1]
        unique_base = split(unique_result[14:end], ".")[1]  # Remove prefix (13 chars + 1)
        
        # Convert hyphens to colons in time portions for DateTime parsing
        micro_datetime = replace(micro_base, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
        nano_datetime = replace(nano_base, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
        precise_datetime = replace(precise_base, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
        unique_datetime = replace(unique_base, r"T(\d{2})-(\d{2})-(\d{2})" => s"T\1:\2:\3")
        
        # All base timestamps should be very similar (within same second)
        # Parse them to compare
        micro_time = DateTime(micro_datetime, "yyyy-mm-ddTHH:MM:SS")
        nano_time = DateTime(nano_datetime, "yyyy-mm-ddTHH:MM:SS")
        precise_time = DateTime(precise_datetime, "yyyy-mm-ddTHH:MM:SS")
        unique_time = DateTime(unique_datetime, "yyyy-mm-ddTHH:MM:SS")
        
        # All should be within the same second (allowing for execution time)
        micro_ms = abs(Millisecond(micro_time - nano_time).value)
        precise_ms = abs(Millisecond(precise_time - micro_time).value) 
        unique_ms = abs(Millisecond(unique_time - micro_time).value)
        @test micro_ms < 2000  # Within 2 seconds (2000 milliseconds)
        @test precise_ms < 2000
        @test unique_ms < 2000
    end
    
    @testset "Edge cases and robustness" begin
        # Test rapid successive calls don't cause errors
        @test_nowarn begin
            for _ in 1:1000
                now_microseconds()
            end
        end
        
        # Test that functions work with different system loads
        # This is more of a smoke test
        results = String[]
        for _ in 1:20
            push!(results, unique_name())
            # Add some computational load
            sum(rand(100))
        end
        
        # All results should still be unique and properly formatted
        @test length(unique(results)) == 20
        @test all(r -> occursin(r"^floridyn_run_\d{4}-\d{2}-\d{2}T\d{2}-\d{2}-\d{2}\.\d{6}$", r), results)
    end
    
    @testset "delete_results function" begin
        @testset "Input validation" begin
            # Create a test Vis object
            vis = Vis(online=false, output_folder="test_validation")
            vis.unique_folder = ""  # Empty for these tests
            
            # Test with non-positive n
            @test_logs (:warn, r"Number of folders to delete must be positive") begin
                result = delete_results(vis, 0)
                @test isempty(result)
            end
            
            @test_logs (:warn, r"Number of folders to delete must be positive") begin
                result = delete_results(vis, -1)
                @test isempty(result)
            end
            
            # Clean up
            rm(joinpath(pwd(), "test_validation"), recursive=true, force=true)
        end
        
        @testset "Empty directory behavior" begin
            # Create a test Vis object
            vis = Vis(online=false, output_folder="test_empty")
            vis.unique_folder = ""  # Empty for these tests
            
            # Test with empty directory (no floridyn_run folders)
            @test_logs (:info, r"No floridyn_run directories found") begin
                result = delete_results(vis, 1)
                @test isempty(result)
            end
            
            # Test dry run with empty directory
            @test_logs (:info, r"No floridyn_run directories found") begin
                result = delete_results(vis, 1, dry_run=true)
                @test isempty(result)
            end
            
            # Clean up
            rm(joinpath(pwd(), "test_empty"), recursive=true, force=true)
        end

        @testset "Single folder operations" begin
            # Create a test Vis object
            vis = Vis(online=false, output_folder="test_single")
            vis.unique_folder = ""  # Empty for these tests
            output_dir = vis.output_path  # This creates the directory structure
            
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
            
            # Test vis.unique_folder is cleared
            @test vis.unique_folder == ""
            
            # Clean up
            rm(joinpath(pwd(), "test_single"), recursive=true, force=true)
        end
            
        @testset "Multiple folder operations" begin
            # Create a test Vis object 
            vis = Vis(online=false, output_folder="test_multiple")
            vis.unique_folder = ""  # Empty for these tests
            output_dir = vis.output_path
            
            # Create multiple test directories with different timestamps
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
                @test length(result) == 2
                # Should return newest first
                @test basename(result[1]) == "floridyn_run_2025-01-01T12-00-00.333333"
                @test basename(result[2]) == "floridyn_run_2025-01-01T11-00-00.222222"
                # All folders should still exist
                @test all(isdir, folder_paths)
            end
            
            # Test actual deletion of 2 newest
            @test_logs (:info, r"Deleting") match_mode=:any begin
                result = delete_results(vis, 2)
                @test length(result) == 2
                # Should delete newest first
                @test basename(result[1]) == "floridyn_run_2025-01-01T12-00-00.333333"
                @test basename(result[2]) == "floridyn_run_2025-01-01T11-00-00.222222"
                
                # Check which folders remain
                @test !isdir(folder_paths[3])  # Newest deleted
                @test !isdir(folder_paths[2])  # Middle deleted
                @test isdir(folder_paths[1])   # Oldest preserved
            end
            
            # Clean up
            rm(joinpath(pwd(), "test_multiple"), recursive=true, force=true)
        end
        
        @testset "Boundary conditions" begin
            # Create a test Vis object
            vis = Vis(online=false, output_folder="test_boundary")
            vis.unique_folder = ""  # Empty for these tests
            output_dir = vis.output_path
            
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
                result = delete_results(vis, 10)  # Request 10, only 3 available
                @test length(result) == 3  # Should only delete what's available
                @test all(!isdir, folder_paths)  # All should be deleted
            end
            
            # Clean up
            rm(joinpath(pwd(), "test_boundary"), recursive=true, force=true)
        end
        
        @testset "Error handling" begin
            # Create a test Vis object
            vis = Vis(online=false, output_folder="test_error")
            vis.unique_folder = ""  # Empty for these tests
            output_dir = vis.output_path
            
            # Create a test directory
            test_folder = joinpath(output_dir, "floridyn_run_2025-01-01T12-00-00.123456")
            mkdir(test_folder)
            
            # Test deletion - should handle errors gracefully
            # Note: This might succeed on some systems, so we mainly test it doesn't crash
            result = delete_results(vis, 1)
            @test isa(result, Vector{String})  # Should return a vector
            
            # Clean up
            rm(joinpath(pwd(), "test_error"), recursive=true, force=true)
        end
        
        @testset "Integration with find_floridyn_runs" begin
            # Create a test Vis object
            vis = Vis(online=false, output_folder="test_integration")
            vis.unique_folder = ""  # Empty for these tests
            output_dir = vis.output_path
            
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
            @test length(result) == 3  # Only the 3 floridyn_run directories
            
            # Check that non-floridyn directories are preserved
            @test isdir(joinpath(output_dir, "other_folder"))
            @test isdir(joinpath(output_dir, "not_floridyn_folder"))
            
            # Check that floridyn_run directories are deleted
            @test !isdir(joinpath(output_dir, "floridyn_run_2025-01-01T10-00-00.111111"))
            @test !isdir(joinpath(output_dir, "floridyn_run_2025-01-01T11-00-00.222222"))
            @test !isdir(joinpath(output_dir, "floridyn_run_2025-01-01T12-00-00.333333"))
            
            # Clean up
            rm(joinpath(pwd(), "test_integration"), recursive=true, force=true)
        end
        
        @testset "Default parameter behavior" begin
            # Create a test Vis object
            vis = Vis(online=false, output_folder="test_default")
            vis.unique_folder = ""  # Empty for these tests
            output_dir = vis.output_path
            
            # Test default n=1
            test_folder = joinpath(output_dir, "floridyn_run_2025-01-01T12-00-00.123456")
            mkdir(test_folder)
            
            # Test that calling without explicit n deletes 1 folder
            result = delete_results(vis)  # Uses default n=1
            @test length(result) == 1
            @test !isdir(test_folder)
            
            # Clean up
            rm(joinpath(pwd(), "test_default"), recursive=true, force=true)
        end
    end
end
