# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
Interactive test runner for FLORIDyn.jl

This script provides a menu-driven interface to run individual test files
or groups of tests from the test/ directory.
"""

# ── Enforce system Python if install_controlplots selected --system ───────────
let prefs_file = joinpath(dirname(@__DIR__), "LocalPreferences.toml")
    if isfile(prefs_file)
        lines = readlines(prefs_file; keep=true)
        in_condapkg = false
        in_pythoncall = false
        backend_null = false
        python_exe = ""
        for line in lines
            if startswith(line, "[CondaPkg]")
                in_condapkg = true
                in_pythoncall = false
            elseif startswith(line, "[PythonCall]")
                in_pythoncall = true
                in_condapkg = false
            elseif startswith(line, "[")
                in_condapkg = false
                in_pythoncall = false
            elseif in_condapkg && contains(line, "backend") && contains(line, "Null")
                backend_null = true
            elseif in_pythoncall && contains(line, "exe")
                python_exe = strip(split(line, '=')[2], [' ', '"', '\t', '\n', '\r'])
            end
        end
        if backend_null
            ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
            if !isempty(python_exe)
                ENV["JULIA_PYTHONCALL_EXE"] = python_exe
            end
        end
    end
end

using Pkg
if ! ("Aqua" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using FLORIDyn
using Test
using LinearAlgebra
using Random
using Suppressor

if basename(pwd()) == "test"
    cd("..")
end

# Files that need error suppression
suppress_error_files = [
    "test_dir.jl",
    "test_shear.jl",
    "test_tit.jl", 
    "test_vel.jl",
    "test_iterate.jl",
    "test_floridyn_cl.jl",
    "test_prepare_simulation.jl"
]

function get_test_files()
    """Get all test_*.jl files from the test directory, excluding helper and empty files"""
    test_dir = @__DIR__
    files = readdir(test_dir)
    
    # Files to exclude from the menu (helper files, empty files, benchmark files)
    excluded_files = ["test_helpers.jl", "test_init.jl", "test_floris2.jl"]
    
    test_files = filter(f -> startswith(f, "test_") && endswith(f, ".jl") && !(f in excluded_files), files)
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
    println("\n" * "="^66)
    println("Running: $description")
    println("File: $filepath")
    println("="^66); println()
    
    # Run the file
    success = false
    try   
        if basename(filepath) in suppress_error_files
            @suppress_err include(filepath)
        else
            include(filepath)
        end
        println("\n✅ $description completed successfully!")
        success = true
    catch e
        println("\n❌ $description failed with error:")
        println(e)
        success = false
    end
    
    return success
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
            if choice in [11, 18]
                  Pkg.activate(dirname(@__DIR__))
            elseif choice == 26
                if ! ("Aqua" ∈ keys(Pkg.project().dependencies))
                    TestEnv.activate()
                end
            end
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
            @info "Running all tests..."
            @eval Main using Pkg
            Main.Pkg.activate(dirname(@__DIR__))
            Main.Pkg.test()
            
        elseif choice == special_start + 1
            # Run Aqua.jl
            filepath = joinpath(@__DIR__, "aqua.jl")
            println(pwd())
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
