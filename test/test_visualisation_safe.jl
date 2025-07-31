# Safe Visualization Testing for Interactive Sessions
# Copyright (c) 2025 Uwe Fechner  
# SPDX-License-Identifier: BSD-3-Clause

"""
    run_visualization_tests()

Safely run visualization tests in interactive sessions with Revise loaded.
This avoids the segmentation fault that occurs when using include("test/test_visualisation.jl").
"""
function run_visualization_tests()
    println("🧪 Running visualization tests safely in interactive session...")
    
    # Ensure Test is available 
    @eval Main using Test
    
    # Run tests via Pkg.test (safest approach)
    @eval Main using Pkg
    Main.Pkg.test(test_args=["test_visualisation.jl"])
end

run_visualization_tests()