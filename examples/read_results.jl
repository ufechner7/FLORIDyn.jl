# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Script to read simulation results from the last run, if it was saved as "results.jld2".
# Creates global variables: wf, md, mi

using JLD2
using FLORIDyn

function read_results(filepath::String = "out/results.jld2")
    # Auto-discover newest run if default path not found and user passed the default
    if filepath == "out/results.jld2" && !isfile(filepath)
        newest = find_newest_run("out")
        if newest !== nothing
            candidate = joinpath("out", newest, "results.jld2")
            if isfile(candidate)
                filepath = candidate
            end
        end
    end
    # Final existence check
    if !isfile(filepath)
        error("Results file not found: $filepath")
    end
    
    @info "Loading simulation results from: $filepath"
    
    try
        # Load the JLD2 file
        data = JLD2.load(filepath)
        
        # Check if required variables are present
        required_vars = ["wf", "md", "mi"]
        missing_vars = String[]
        
        for var in required_vars
            if !haskey(data, var)
                push!(missing_vars, var)
            end
        end
        
        if !isempty(missing_vars)
            error("Missing required variables in results file: $(join(missing_vars, ", "))")
        end
        
        # Extract the variables
        wf = data["wf"]
        md = data["md"]
        mi = data["mi"]
                
        @info "Successfully loaded simulation results:"
        @info "  - WindFarm (wf): $(typeof(wf))"
        @info "  - Measurement data   (md): $(typeof(md))"
        @info "  - Interaction matrix (mi): $(typeof(mi))"
        
        return wf, md, mi
        
    catch e
        if isa(e, SystemError)
            error("Failed to read file $filepath: $(e.msg)")
        else
            error("Failed to load simulation results: $e")
        end
    end
end

vis = Vis("data/vis_default.yaml")
runs = find_floridyn_runs(vis.output_folder)
if length(runs) == 0
    @warn "No FLORIDyn runs found in output folder: $(vis.output_folder)"
else
    wf, md, mi = read_results(joinpath(runs[end], "results.jld2"))
end

nothing