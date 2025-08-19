#!/usr/bin/env julia
# Determine minimal and maximal x and y positions of wind turbines from a YAML config.
# Usage:
#   julia scripts/minmax_turbine_positions.jl [path/to/config.yaml]
# If no path is provided, defaults to data/2021_54T_NordseeOne.yaml relative to repo root.

using YAML

function minmax_turbine_positions(yaml_path::String)
    cfg = YAML.load_file(yaml_path)
    haskey(cfg, "turbines") || error("No 'turbines' key found in: $(yaml_path)")
    tlist = cfg["turbines"]
    isempty(tlist) && error("Empty turbines list in: $(yaml_path)")

    xs = Float64[]
    ys = Float64[]
    ids = Int[]
    for t in tlist
        # YAML.jl returns Dict{String,Any} entries
        haskey(t, "x") || error("Turbine entry missing 'x': $(t)")
        haskey(t, "y") || error("Turbine entry missing 'y': $(t)")
        push!(xs, float(t["x"]))
        push!(ys, float(t["y"]))
        push!(ids, Int(t["id"]))
    end
    minx, maxx = minimum(xs), maximum(xs)
    miny, maxy = minimum(ys), maximum(ys)
    (; minx, maxx, miny, maxy, n=length(xs), ids_extrema=(ids[argmin(xs)], ids[argmax(xs)], ids[argmin(ys)], ids[argmax(ys)]))
end

function print_report(res, path)
    println("File: ", path)
    println("Turbines: ", res.n)
    println("x_min: ", res.minx)
    println("x_max: ", res.maxx)
    println("y_min: ", res.miny)
    println("y_max: ", res.maxy)
    xmin_id, xmax_id, ymin_id, ymax_id = res.ids_extrema
    println("x_min_id: ", xmin_id, ", x_max_id: ", xmax_id)
    println("y_min_id: ", ymin_id, ", y_max_id: ", ymax_id)
end

function main()
    default_path = joinpath(@__DIR__, "..", "data", "2021_54T_NordseeOne.yaml")
    path = get(ARGS, 1, default_path)
    res = minmax_turbine_positions(abspath(path))
    print_report(res, abspath(path))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
