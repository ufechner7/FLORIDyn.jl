#!/usr/bin/env julia
# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
Interactive menu for FLORIDyn.jl examples and documentation.

This script provides a simple menu interface to run examples and access documentation
for the FLORIDyn.jl wind farm simulation package.
"""

using Pkg
using FLORIDyn
using REPL.TerminalMenus

settings_file, vis_file = get_default_project()
@info "Using settings file: $settings_file"

options = ["\"select_project()\"; select_project()",
           "\"flow_field_vel_reduction\"; PLT=1; include(\"main.jl\")",
           "\"flow_field_added_turbulence\"; PLT=2; include(\"main.jl\")",
           "\"flow_field_eff_wind_speed\"; PLT=3; include(\"main.jl\")",
           "\"plot_measurements\"; PLT=4; include(\"main.jl\")",
           "\"plot_measurements_lineplot\"; PLT=5; include(\"main.jl\")",
           "\"flow_field_vel_reduction_online\"; PLT=6; include(\"main.jl\")",
           "\"create_video_from_saved_frames\"; PLT=7; include(\"main.jl\")",
           "\"run_all_visualisations\"; include(\"main_all.jl\")",
           "\"read_results\"; include(\"read_results.jl\")",
           "\"plot_wind_direction\"; include(\"plot_wind_dir.jl\")",
           "\"play_videos\"; include(\"play_video.jl\")",
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
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

# Main entry point
example_menu()
