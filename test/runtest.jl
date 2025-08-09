#!/usr/bin/env julia

"""
Interactive test runner for FLORIDyn.jl

This script provides a menu-driven interface to run individual test files
or groups of tests from the test/ directory.
"""

using Pkg
if ! ("Suppressor" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using FLORIDyn
using Test
using LinearAlgebra
using Random
using Suppressor
using DistributedNext

if basename(pwd()) == "test"
    cd("..")
end

function get_test_files()
    """Get all test_*.jl files from the test directory"""
    test_dir = @__DIR__
    files = readdir(test_dir)
    test_files = filter(f -> startswith(f, "test_") && endswith(f, ".jl"), files)
    return sort(test_files)
end

function get_bench_files()
    """Get all bench_*.jl files from the test directory"""
    test_dir = @__DIR__
    files = readdir(test_dir)
    bench_files = filter(f -> startswith(f, "bench_") && endswith(f, ".jl"), files)
    return sort(bench_files)
end

function display_menu()
    """Display the main menu"""
    println("\n" * "="^60)
    println("         FLORIDyn.jl Test Runner")
    println("="^60)
    
    # Get test files
    test_files = get_test_files()
    bench_files = get_bench_files()
    
    println("\nSelect a test to run:")
    println()
    
    # Display test files
    println("📋 Test Files:")
    for (i, file) in enumerate(test_files)
        test_name = replace(file, "test_" => "", ".jl" => "")
        println("  $i. $test_name")
    end
    
    println()
    
    # Display benchmark files
    if !isempty(bench_files)
        println("⚡ Benchmark Files:")
        start_idx = length(test_files) + 1
        for (i, file) in enumerate(bench_files)
            bench_name = replace(file, "bench_" => "", ".jl" => "")
            println("  $(start_idx + i - 1). $bench_name (benchmark)")
        end
        println()
    end
    
    # Special options
    special_start = length(test_files) + length(bench_files) + 1
    println("📦 Special Options:")
    println("  $(special_start). Run all tests (runtests.jl)")
    println("  $(special_start + 1). Run Aqua.jl quality checks")
    println("  $(special_start + 2). Exit")
    println()
    
    return test_files, bench_files, special_start
end

function run_file(filepath::String, description::String)
    """Run a specific test file"""
    println("\n" * "="^60)
    println("Running: $description")
    println("File: $filepath")
    println("="^60)
    
    # Activate the project environment from the parent directory (project root)
    Pkg.activate(dirname(@__DIR__))
    
    # Run the file
    try
        include(filepath)
        println("\n✅ $description completed successfully!")
    catch e
        println("\n❌ $description failed with error:")
        println(e)
        return false
    end
    return true
end

function main()
    """Main interactive loop"""
    while true
        test_files, bench_files, special_start = display_menu()
        
        print("Enter your choice (1-$(special_start + 2)): ")
        input = readline()
        
        # Parse input
        local choice
        try
            choice = parse(Int, strip(input))
        catch
            println("❌ Invalid input. Please enter a number.")
            continue
        end
        
        # Handle choice
        total_tests = length(test_files)
        total_bench = length(bench_files)
        
        if choice >= 1 && choice <= total_tests
            # Regular test file
            file = test_files[choice]
            test_name = replace(file, "test_" => "", ".jl" => "")
            filepath = joinpath(@__DIR__, file)
            run_file(filepath, "Test: $test_name")
            
        elseif choice >= total_tests + 1 && choice <= total_tests + total_bench
            # Benchmark file
            bench_idx = choice - total_tests
            file = bench_files[bench_idx]
            bench_name = replace(file, "bench_" => "", ".jl" => "")
            filepath = joinpath(@__DIR__, file)
            run_file(filepath, "Benchmark: $bench_name")
            
        elseif choice == special_start
            # Run all tests
            filepath = joinpath(@__DIR__, "runtests.jl")
            run_file(filepath, "All Tests")
            
        elseif choice == special_start + 1
            # Run Aqua.jl
            filepath = joinpath(@__DIR__, "aqua.jl")
            run_file(filepath, "Aqua.jl Quality Checks")
            
        elseif choice == special_start + 2
            # Exit
            println("👋 Goodbye!")
            break
            
        else
            println("❌ Invalid choice. Please select a number from the menu.")
            continue
        end
        
        # Ask if user wants to continue
        println("\nPress Enter to return to menu, or 'q' to quit...")
        input = readline()
        if lowercase(strip(input)) == "q"
            println("👋 Goodbye!")
            break
        end
    end
end

# Run the main function when this script is included
main()
