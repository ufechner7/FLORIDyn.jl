# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Script to read simulation results from results.jld2 file
# Creates global variables: wf, md, mi

using JLD2
using FLORIDyn

"""
    find_newest_run(base::String = "out") -> Union{String, Nothing}

Return the name (relative to `base`) of the newest directory whose name starts with
"floridyn_run" inside `base`. If no matching directory exists, return `nothing`.

Selection logic:
1. List entries in `base` that are directories and start with "floridyn_run".
2. Sort them by modification time (mtime) descending.
3. Return the first (newest) name.

Edge cases handled:
- Non‑existent base directory -> returns `nothing`.
- No matching directories -> returns `nothing`.
- Permission errors are propagated as an error.
"""
function find_newest_run(base::String = "out")::Union{String, Nothing}
    isdir(base) || return nothing
    entries = readdir(base; join=false)
    runs = filter(name -> startswith(name, "floridyn_run") && isdir(joinpath(base, name)), entries)
    isempty(runs) && return nothing
    sort!(runs; by = name -> stat(joinpath(base, name)).mtime, rev=true)
    return first(runs)
end

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
        data = load(filepath)
        
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

try
    wf, md, mi = read_results()
catch e
    @warn "Could not auto-load results: $(e)" 
end
nothing