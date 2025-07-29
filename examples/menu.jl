#!/usr/bin/env julia
# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
Interactive menu for FLORIDyn.jl examples and documentation.

This script provides a simple menu interface to run examples and access documentation
for the FLORIDyn.jl wind farm simulation package.
"""

using Pkg
using REPL.TerminalMenus

options = ["main_example = include(\"main.jl\")",
           "open_documentation()",
           "quit"]

function open_documentation()
    println("\nOpening documentation in browser...")
    doc_url = "https://ufechner7.github.io/FLORIDyn.jl/dev/"
    
    try
        if Sys.islinux()
            run(pipeline(`xdg-open $doc_url`, stdout=devnull, stderr=devnull))
        elseif Sys.iswindows()
            run(pipeline(`cmd /c start $doc_url`, stdout=devnull, stderr=devnull))
        elseif Sys.isapple()
            run(pipeline(`open $doc_url`, stdout=devnull, stderr=devnull))
        else
            println("Cannot automatically open browser on this system.")
            println("Please manually open: $doc_url")
        end
        println("Documentation URL: $doc_url")
    catch e
        println("Could not open browser automatically.")
        println("Please manually open: $doc_url")
    end
end

function example_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            try
                eval(Meta.parse(options[choice]))
            catch e
                println("Error executing: ", options[choice])
                println("Error: $e")
                if choice == 1
                    println("Make sure you are in the FLORIDyn.jl root directory and the examples/main.jl file exists.")
                end
            end
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

# Main entry point
example_menu()
