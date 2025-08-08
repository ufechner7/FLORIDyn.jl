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
end
