#!/usr/bin/env julia
# Extend SOWFA_nacelleYaw.csv blocks from 9 turbines (0..8) to 54 turbines (0..53)
# Usage: julia scripts/expand_nacelle_yaw.jl [input_csv] [output_csv]
# If output_csv is omitted, the input file is overwritten in-place.

function expand_nacelle_yaw(input_path::String, output_path::String=input_path)
    lines = readlines(input_path)
    isempty(lines) && error("Empty file: $(input_path)")

    header = lines[1]
    blocks = Vector{Vector{String}}()
    current = String[]

    for ln in lines[2:end]
        s = strip(ln)
        if isempty(s)
            if !isempty(current)
                push!(blocks, current)
                current = String[]
            else
                # multiple blank lines, keep one
                push!(blocks, String[])
            end
        else
            push!(current, ln)
        end
    end
    if !isempty(current)
        push!(blocks, current)
    end

    out = IOBuffer()
    println(out, header)

    n_blocks_written = 0
    for blk in blocks
        if isempty(blk)
            println(out) # preserve empty separator
            continue
        end
        # Build a map tid->line, also validate time/dt consistency and capture reference yaw
        tid2line = Dict{Int, String}()
        ref_time = nothing
        ref_dt = nothing
        ref_yaw = nothing
        for ln in blk
            parts = split(strip(ln))
            length(parts) >= 4 || error("Unexpected line format: '$(ln)'")
            tid = try
                parse(Int, parts[1])
            catch
                error("Unexpected turbine id in line: '$(ln)'")
            end
            time_s = parts[2]
            dt_s   = parts[3]
            yaw_s  = parts[4]
            if ref_time === nothing
                ref_time = time_s
                ref_dt   = dt_s
                ref_yaw  = yaw_s
            else
                if (time_s != ref_time) || (dt_s != ref_dt)
                    error("Inconsistent time/dt within a block: '$(blk[1])' vs '$(ln)'")
                end
            end
            tid2line[tid] = ln
        end
        # Prefer yaw of turbine 0 if available
        if haskey(tid2line, 0)
            parts0 = split(strip(tid2line[0]))
            ref_yaw = parts0[4]
        end

        # Emit turbines in order 0..53; keep existing lines exactly; synthesize missing
        for tid in 0:53
            if haskey(tid2line, tid)
                println(out, tid2line[tid])
            else
                println(out, string(tid, " ", ref_time, " ", ref_dt, " ", ref_yaw))
            end
        end
        println(out) # blank line after each block
        n_blocks_written += 1
    end

    # If overwriting, create a simple backup alongside first
    if output_path == input_path
        backup = input_path * ".bak"
        cp(input_path, backup; force=true)
    end

    open(output_path, "w") do io
        write(io, String(take!(out)))
    end
    return n_blocks_written
end

function main()
    input_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "data", "2021_54T_NordseeOne", "SOWFA_nacelleYaw.csv")
    output_path = length(ARGS) >= 2 ? ARGS[2] : input_path
    nb = expand_nacelle_yaw(abspath(input_path), abspath(output_path))
    println("Wrote ", nb, " blocks to ", output_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
