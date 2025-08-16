# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

module TestHelpers
using FLORIDyn, DataFrames, MAT

export wf_dict2windfarm, structs_equal, compare_windFarms

"""
    wf_dict2windfarm(wf_dict::Dict{String, Any}) -> WindFarm

Convert MATLAB wf_dict to Julia WindFarm type.

This function handles the conversion of WindFarm data structures from MATLAB MAT files
to the Julia WindFarm type. It properly handles the nested cell arrays, mixed data types,
and different array orientations that result from MATLAB to Julia conversion.

# Arguments
- `wf_dict::Dict{String, Any}`: Dictionary from MAT file containing WindFarm data with keys:
  - `"nT"`: Number of turbines (scalar)
  - `"nOP"`: Number of operating points (scalar) 
  - `"States_WF"`: Wind field states matrix
  - `"States_OP"`: Operating point states matrix
  - `"States_T"`: Turbine states matrix
  - `"posBase"`: Base positions matrix (nT × 3)
  - `"posNac"`: Nacelle positions matrix (nT × 3) 
  - `"D"`: Turbine diameters (nT × 1)
  - `"StartI"`: Start indices (1 × nT)
  - `"intOPs"`: Interpolated operating points (nT × 1 cell array)
  - `"weight"`: Weights (nT × 1 cell array, may contain scalars or matrices)
  - `"dep"`: Dependencies (nT × 1 cell array, may contain scalars, vectors or matrices)
  - `"red_arr"`: Reduced array matrix
  - `"Names_T"`: Turbine state names (1 × n_states)
  - `"Names_WF"`: Wind field names (1 × n_wf)
  - `"Names_OP"`: Operating point names (1 × n_op)

# Returns
- `WindFarm`: Properly constructed WindFarm object with all fields correctly typed

# Examples
```julia
# Load MAT file and convert
vars = matread("wind_farm_data.mat")
wf_dict = vars["wf"]
wf = wf_dict2windfarm(wf_dict)

# Access as DataFrames
turbine_data = wf.turbines
windfield_data = wf.windfield  
ops_data = wf.ops
```

# Notes
- Handles mixed scalar/array types in weight and dep fields
- Transposes position matrices to match Julia convention (3 × nT)
- Converts MATLAB row vectors to Julia column vectors for names
- Ensures all nested structures are properly typed for WindFarm constructor
"""
function wf_dict2windfarm(wf_dict::Dict{String, Any})
    # Extract and convert the complex nested structures
    nT = Int(wf_dict["nT"])
    nOP = Int(wf_dict["nOP"])
    
    # Convert nested cell arrays
    intOPs_vec = Vector{Matrix{Float64}}()
    weight_vec = Vector{Vector{Float64}}()
    dep_vec = Vector{Vector{Int}}()
    
    for i in 1:nT
        if  "intOPs" in keys(wf_dict)
            push!(intOPs_vec, wf_dict["intOPs"][i])
        end
        
        # Convert weight (may be empty matrix, scalar, or vector)
        if  "weight" in keys(wf_dict)
            w = wf_dict["weight"][i]
            if isempty(w)
                push!(weight_vec, Float64[])
            elseif isa(w, Number)  # Handle scalar case
                push!(weight_vec, [Float64(w)])
            else
                push!(weight_vec, vec(Float64.(w)))
            end
        end
        if  "dep" in keys(wf_dict)
            # Convert dep (may be empty matrix, scalar, or vector)
            d = wf_dict["dep"][i]
            if isempty(d)
                push!(dep_vec, Int[])
            elseif isa(d, Number)  # Handle scalar case
                push!(dep_vec, [Int(d)])
            else
                push!(dep_vec, vec(Int.(d)))
            end
        end
    end
    
    return WindFarm(
        nT = nT,
        nOP = nOP,
        States_WF = wf_dict["States_WF"],
        States_OP = wf_dict["States_OP"], 
        States_T = wf_dict["States_T"],
        posBase = wf_dict["posBase"],
        posNac = wf_dict["posNac"],
        D = vec(wf_dict["D"]),
        StartI = Int.(wf_dict["StartI"]),
        intOPs = intOPs_vec,
        Weight = weight_vec,
        dep = dep_vec,
        red_arr = wf_dict["red_arr"],
        Names_T = String.(wf_dict["Names_T"][1, :]),      # Get row vector and convert to strings
        Names_WF = String.(wf_dict["Names_WF"][1, :]),    # Get row vector and convert to strings
        Names_OP = String.(wf_dict["Names_OP"][1, :])     # Get row vector and convert to strings
    )
end

function structs_equal(a::T, b::T; prn=true) where T
    result = true
    fields = fieldnames(T)
    for f in fields
        val_a = getfield(a, f)
        val_b = getfield(b, f)
        if val_a != val_b
            prn && println("Field $(f): a = $(val_a), b = $(val_b)")
            result = false
        end
    end
    return result
end

"""
    compare_windFarms(wf1::WindFarm, wf2::WindFarm; detailed=true, tolerance=1e-10) -> Bool

Compare two WindFarm objects and print detailed differences between them.

This function provides a comprehensive comparison of WindFarm objects, showing
differences in dimensions, field values, and data content with optional tolerance
for floating-point comparisons.

# Arguments
- `wf1::WindFarm`: First WindFarm object to compare
- `wf2::WindFarm`: Second WindFarm object to compare  
- `detailed::Bool=true`: Whether to show detailed field-by-field comparison
- `tolerance::Float64=1e-10`: Tolerance for floating-point comparisons

# Returns
- `Bool`: `true` if WindFarms are equal (within tolerance), `false` otherwise

# Examples
```julia
wf1 = WindFarm(nT=3, nOP=100, ...)
wf2 = WindFarm(nT=3, nOP=100, ...)

# Basic comparison
are_equal = compare_windFarms(wf1, wf2)

# Comparison with custom tolerance  
are_equal = compare_windFarms(wf1, wf2, tolerance=1e-8)

# Silent comparison
are_equal = compare_windFarms(wf1, wf2, detailed=false)
```
"""
function compare_windFarms(wf1::WindFarm, wf2::WindFarm; detailed=true, tolerance=1e-10)
    if detailed
        println("=" ^ 60)
        println("WindFarm Comparison")
        println("=" ^ 60)
    end
    
    all_equal = true
    
    # Compare basic dimensions
    if detailed
        println("\n📊 Basic Dimensions:")
        println("  Field                WF1        WF2        Equal")
        println("  " * "-" ^ 45)
    end
    
    basic_fields = [:nT, :nOP]
    for field in basic_fields
        val1 = getfield(wf1, field)
        val2 = getfield(wf2, field)
        equal = (val1 == val2)
        all_equal &= equal
        
        if detailed
            status = equal ? "✓" : "✗"
            println(@sprintf("  %-15s %10s %10s   %s", string(field), val1, val2, status))
        end
    end
    
    # Compare matrix/vector dimensions
    if detailed
        println("\n📏 Matrix/Vector Dimensions:")
        println("  Field                WF1 Size      WF2 Size      Equal")
        println("  " * "-" ^ 55)
    end
    
    matrix_fields = [:States_WF, :States_OP, :States_T, :posBase, :posNac, :D, :StartI, :red_arr]
    for field in matrix_fields
        val1 = getfield(wf1, field)
        val2 = getfield(wf2, field)
        size1 = size(val1)
        size2 = size(val2)
        equal = (size1 == size2)
        all_equal &= equal
        
        if detailed
            status = equal ? "✓" : "✗"
            println(sprintf("  %-15s %-12s %-12s   %s", string(field), string(size1), string(size2), status))
        end
    end
    
    # Compare vector lengths for nested structures
    if detailed
        println("\n📋 Vector Lengths:")
        println("  Field                WF1 Length   WF2 Length   Equal")
        println("  " * "-" ^ 50)
    end
    
    vector_fields = [:Names_T, :Names_WF, :Names_OP, :intOPs, :Weight, :dep]
    for field in vector_fields
        val1 = getfield(wf1, field)
        val2 = getfield(wf2, field)
        len1 = length(val1)
        len2 = length(val2)
        equal = (len1 == len2)
        all_equal &= equal
        
        if detailed
            status = equal ? "✓" : "✗"
            println(sprintf("  %-15s %10d %10d     %s", string(field), len1, len2, status))
        end
    end
    
    # Compare string vectors (names)
    if detailed
        println("\n🏷️  Name Vectors:")
    end
    
    name_fields = [:Names_T, :Names_WF, :Names_OP]
    for field in name_fields
        val1 = getfield(wf1, field)
        val2 = getfield(wf2, field)
        equal = (val1 == val2)
        all_equal &= equal
        
        if detailed && !equal
            println("  $(field) differs:")
            println("    WF1: $(val1)")
            println("    WF2: $(val2)")
        elseif detailed
            println("  $(field): ✓ ($(length(val1)) elements)")
        end
    end
    
    # Compare numerical matrices with tolerance
    if detailed
        println("\n🔢 Numerical Data (tolerance=$(tolerance)):")
    end
    
    numerical_fields = [:States_WF, :States_OP, :States_T, :posBase, :posNac, :D, :StartI, :red_arr]
    for field in numerical_fields
        val1 = getfield(wf1, field)
        val2 = getfield(wf2, field)
        
        if size(val1) == size(val2)
            if isempty(val1) && isempty(val2)
                equal = true
            else
                max_diff = maximum(abs.(val1 .- val2))
                equal = max_diff <= tolerance
            end
            all_equal &= equal
            
            if detailed
                if equal
                    println("  $(field): ✓")
                else
                    println("  $(field): ✗ (max difference: $(max_diff))")
                end
            end
        end
    end
    
    # Compare nested structures  
    if detailed
        println("\n🎯 Nested Structures:")
    end
    
    # Compare intOPs
    intops_equal = true
    if length(wf1.intOPs) == length(wf2.intOPs)
        for i in 1:length(wf1.intOPs)
            if size(wf1.intOPs[i]) != size(wf2.intOPs[i])
                intops_equal = false
                break
            elseif !isempty(wf1.intOPs[i]) && maximum(abs.(wf1.intOPs[i] .- wf2.intOPs[i])) > tolerance
                intops_equal = false
                break
            end
        end
    else
        intops_equal = false
    end
    all_equal &= intops_equal
    
    if detailed
        status = intops_equal ? "✓" : "✗"
        println("  intOPs: $(status)")
    end
    
    # Compare Weight
    weight_equal = true
    if length(wf1.Weight) == length(wf2.Weight)
        for i in 1:length(wf1.Weight)
            if length(wf1.Weight[i]) != length(wf2.Weight[i])
                weight_equal = false
                break
            elseif !isempty(wf1.Weight[i]) && maximum(abs.(wf1.Weight[i] .- wf2.Weight[i])) > tolerance
                weight_equal = false
                break
            end
        end
    else
        weight_equal = false
    end
    all_equal &= weight_equal
    
    if detailed
        status = weight_equal ? "✓" : "✗"
        println("  Weight: $(status)")
    end
    
    # Compare dep
    dep_equal = (wf1.dep == wf2.dep)
    all_equal &= dep_equal
    
    if detailed
        status = dep_equal ? "✓" : "✗"
        println("  dep: $(status)")
    end
    
    # Final summary
    if detailed
        println("\n" * "=" ^ 60)
        if all_equal
            println("🎉 WindFarms are EQUAL")
        else
            println("⚠️  WindFarms are DIFFERENT")
        end
        println("=" ^ 60)
    end
    
    return all_equal
end

# Helper function for formatting (simple sprintf replacement)
function sprintf(fmt, args...)
    # Simple implementation for the format strings we use
    if fmt == "  %-15s %10s %10s   %s"
        return "  $(rpad(args[1], 15)) $(lpad(string(args[2]), 10)) $(lpad(string(args[3]), 10))   $(args[4])"
    elseif fmt == "  %-15s %-12s %-12s   %s"
        return "  $(rpad(args[1], 15)) $(rpad(args[2], 12)) $(rpad(args[3], 12))   $(args[4])"
    elseif fmt == "  %-15s %10d %10d     %s"
        return "  $(rpad(args[1], 15)) $(lpad(string(args[2]), 10)) $(lpad(string(args[3]), 10))     $(args[4])"
    else
        return string(args...)
    end
end
end