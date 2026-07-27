# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Compare the sampling coordinates from data/velocity_sampling.txt with the x positions
# of the wind turbines in data/windFarmProperties.yaml.
# Usage:
#   julia scripts/check_velocity_sampling.jl [path/to/velocity_sampling.txt] [path/to/windFarmProperties.yaml]

using YAML

"""
    read_sampling_coordinates(path)

Read the `coordinates` entry of a TOSCA velocity sampling file and return it as a
`Vector{Float64}`. Also returns `surface_number`, the number of sampling surfaces
declared in the same file (`nothing` if the key is missing).
"""
function read_sampling_coordinates(path::String)
    coordinates = Float64[]
    surface_number = nothing
    for line in eachline(path)
        tokens = split(strip(line))
        isempty(tokens) && continue
        if tokens[1] == "coordinates"
            isempty(coordinates) || error("Duplicate 'coordinates' entry in: $(path)")
            coordinates = parse.(Float64, tokens[2:end])
        elseif tokens[1] == "surfaceNumber"
            surface_number = parse(Int, tokens[2])
        end
    end
    isempty(coordinates) && error("No 'coordinates' entry found in: $(path)")
    (; coordinates, surface_number)
end

"""
    read_turbine_x(path)

Read the turbine x positions from a wind farm properties YAML file. Returns the
turbine names (sorted naturally by turbine number) and the corresponding x positions.
"""
function read_turbine_x(path::String)
    cfg = YAML.load_file(path)
    haskey(cfg, "turbines") || error("No 'turbines' key found in: $(path)")
    turbines = cfg["turbines"]
    isempty(turbines) && error("Empty 'turbines' entry in: $(path)")

    names = sort(collect(keys(turbines)); by = t -> parse(Int, t[2:end]))
    x = [float(turbines[name]["baseLocation"][1]) for name in names]
    (; names, x)
end

function check_velocity_sampling(sampling_path::String, yaml_path::String; atol=1e-6)
    sampling = read_sampling_coordinates(sampling_path)
    farm = read_turbine_x(yaml_path)
    coordinates = sampling.coordinates
    x_turbines = farm.x

    println("Sampling file: ", sampling_path)
    println("Wind farm file: ", yaml_path)
    println("Sampling coordinates: ", length(coordinates))
    println("Turbines: ", length(x_turbines))

    if sampling.surface_number !== nothing && sampling.surface_number != length(coordinates)
        println("WARNING: surfaceNumber ($(sampling.surface_number)) != number of ",
                "coordinates ($(length(coordinates)))")
    end

    if length(coordinates) != length(x_turbines)
        println("MISMATCH: cannot compare element-wise, the counts differ.")
        return (; match=false, offsets=Float64[])
    end

    offsets = coordinates .- x_turbines
    println("\nOffset (sampling coordinate - turbine x):")
    println("  min:  ", minimum(offsets), " m")
    println("  max:  ", maximum(offsets), " m")
    println("  mean: ", sum(offsets) / length(offsets), " m")

    constant = maximum(offsets) - minimum(offsets) <= atol
    identical = all(abs.(offsets) .<= atol)

    if identical
        println("\nOK: all sampling coordinates equal the turbine x positions.")
    elseif constant
        println("\nOK: all sampling coordinates are shifted by a constant offset of ",
                offsets[1], " m relative to the turbine x positions.")
    else
        println("\nMISMATCH: the offset is not constant. Deviations from the mean offset:")
        mean_offset = sum(offsets) / length(offsets)
        for i in eachindex(offsets)
            deviation = offsets[i] - mean_offset
            abs(deviation) <= atol && continue
            println("  ", rpad(farm.names[i], 6), " x = ", rpad(x_turbines[i], 10),
                    " coordinate = ", rpad(coordinates[i], 10),
                    " offset = ", rpad(offsets[i], 10),
                    " deviation = ", deviation)
        end
    end

    (; match=identical || constant, offsets)
end

function main()
    default_sampling = joinpath(@__DIR__, "..", "data", "velocity_sampling.txt")
    default_yaml = joinpath(@__DIR__, "..", "data", "windFarmProperties.yaml")
    sampling_path = abspath(get(ARGS, 1, default_sampling))
    yaml_path = abspath(get(ARGS, 2, default_yaml))
    check_velocity_sampling(sampling_path, yaml_path)
end

main()
