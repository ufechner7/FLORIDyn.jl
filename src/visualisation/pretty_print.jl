# Pretty-printing functions for FLORIDyn structs
# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Printf

"""
    Base.show(io::IO, wf::WindFarm)

Pretty-print a WindFarm struct with formatted output showing key information about the wind farm layout and configuration.

# Example
```julia
julia> show(stdout, wf)
WindFarm:
  Turbines:       9
  Operating Points: 200
  Layout:         3×3 grid (estimated)
  State Names:    a, yaw, TI
  Diameters:      126.0 m (uniform)
  Dependencies:   Calculated for 9 turbines
  
  State Dimensions:
    Wind Field:    7200 states
    Turbines:      5400 states
    Op. Points:   10800 states
  State Names:
    WF: wind_vel, wind_dir, TI0, OP_ori
    OP: x0, y0, z0, x1, y1, z1
```
"""
function Base.show(io::IO, wf::WindFarm)
    println(io, "WindFarm:")
    
    # Basic information
    println(io, "  Turbines:       $(wf.nT)")
    println(io, "  Operating Points: $(wf.nOP)")
    
    # Layout estimation (if positions are available)
    if !isempty(wf.posBase) && size(wf.posBase, 1) > 0
        layout_info = estimate_layout(wf.posBase)
        println(io, "  Layout:         $layout_info")
        
        # Position range
        if size(wf.posBase, 1) >= 2
            x_range = maximum(wf.posBase[1, :]) - minimum(wf.posBase[1, :])
            y_range = maximum(wf.posBase[2, :]) - minimum(wf.posBase[2, :])
            println(io, "  Layout Size:    $(round(x_range, digits=1))m × $(round(y_range, digits=1))m")
        end
    end
      
    # Diameters information
    if !isempty(wf.D)
        if all(d -> d ≈ wf.D[1], wf.D)
            println(io, "  Diameters:      $(wf.D[1]) m (uniform)")
        else
            min_d = minimum(wf.D)
            max_d = maximum(wf.D)
            println(io, "  Diameters:      $(min_d) - $(max_d) m (variable)")
        end
    end
    
    # Dependencies
    if !isempty(wf.dep)
        total_deps = sum(length, wf.dep)
        println(io, "  Dependencies:   $total_deps connections for $(length(wf.dep)) turbines")
    end
    
    # State dimensions
    println(io, "")
    println(io, "  State Dimensions:")
    if !isempty(wf.States_WF)
        println(io, "    Wind Field:    $(size(wf.States_WF)[1]) x $(size(wf.States_WF)[2]) states")
    end
    if !isempty(wf.States_T)
        println(io, "    Turbines:      $(size(wf.States_T)[1]) x $(size(wf.States_T)[2]) states")
    end
    if !isempty(wf.States_OP)
        println(io, "    Op. Points:    $(size(wf.States_OP)[1]) x $(size(wf.States_OP)[2]) states")
    end
    
    # Additional information for very detailed view
    if get(io, :compact, false) == false
        println(io, "")
        println(io, "  State Names:")
        if !isempty(wf.Names_WF)
            if length(wf.Names_WF) <= 10
                println(io, "    WF: $(join(wf.Names_WF, ", "))")
            else
                println(io, "    WF: $(join(wf.Names_WF[1:5], ", "))... (+$(length(wf.Names_WF)-5) more)")
            end
        end
        if !isempty(wf.Names_T)
            if length(wf.Names_T) <= 5
                names_str = join(wf.Names_T, ", ")
            else
                names_str = join(wf.Names_T[1:3], ", ") * ", ..., " * wf.Names_T[end]
            end
            println(io, "     T: $names_str")
        end
        if !isempty(wf.Names_OP)
            if length(wf.Names_OP) <= 10
                println(io, "    OP: $(join(wf.Names_OP, ", "))")
            else
                println(io, "    OP: $(join(wf.Names_OP[1:5], ", "))... (+$(length(wf.Names_OP)-5) more)")
            end
        end
    end
end

"""
    estimate_layout(positions::Matrix{Float64})

Estimate the wind farm layout pattern from turbine positions.
Returns a string description of the estimated layout.
"""
function estimate_layout(positions::Matrix{Float64})
    if size(positions, 2) == 0
        return "Empty"
    end
    
    n_turbines = size(positions, 1)
    
    if n_turbines == 1
        return "Single turbine"
    end
    
    # Try to detect grid patterns
    if size(positions, 1) >= 2
        x_coords = sort(unique(round.(positions[1, :], digits=1)))
        y_coords = sort(unique(round.(positions[2, :], digits=1)))
        
        n_x = length(x_coords)
        n_y = length(y_coords)
        
        # Check if it's a regular grid
        if n_x * n_y == n_turbines
            return "$(n_x)×$(n_y) grid"
        elseif n_x == n_turbines || n_y == n_turbines
            return "Linear arrangement"
        else
            return "Irregular layout"
        end
    else
        return "Unknown layout"
    end
end

"""
    Base.show(io::IO, ::MIME"text/plain", wf::WindFarm)

Detailed pretty-print for WindFarm when displayed in REPL or notebooks.
"""
function Base.show(io::IO, ::MIME"text/plain", wf::WindFarm)
    show(IOContext(io, :compact => false), wf)
end

"""
    summary(wf::WindFarm)

Provide a one-line summary of a WindFarm struct.

# Example
```julia
julia> summary(wf)
"WindFarm: 9 turbines, 3456 op. points, 3×3 grid"
```
"""
function Base.summary(wf::WindFarm)
    layout = ""
    if !isempty(wf.posBase) && size(wf.posBase, 1) > 0
        layout = ", " * estimate_layout(wf.posBase)
    end
    return "WindFarm: $(wf.nT) turbines, $(wf.nOP) op. points$layout"
end

function turbines(wf::WindFarm)
    # Check if we have turbine data
    if isempty(wf.States_T) || isempty(wf.Names_T)
        throw(ArgumentError("WindFarm must have non-empty States_T and Names_T to create DataFrame"))
    end
    
    # Check dimensions match
    if size(wf.States_T, 2) != length(wf.Names_T)
        throw(DimensionMismatch("Number of turbines in States_T ($(size(wf.States_T, 2))) must match length of Names_T ($(length(wf.Names_T)))"))
    end
    
    # Create DataFrame with turbine names as columns
    df = DataFrame()
    
    # Add each turbine's states as a column
    for (i, turbine_name) in enumerate(wf.Names_T)
        df[!, Symbol(turbine_name)] = wf.States_T[:, i]
    end
    
    return df
end

function windfield(wf::WindFarm)
    # Check if we have wind field data
    if isempty(wf.States_WF) || isempty(wf.Names_WF)
        throw(ArgumentError("WindFarm must have non-empty States_WF and Names_WF to create DataFrame"))
    end
    
    # Check dimensions match - States_WF should have as many rows as Names_WF has elements
    if size(wf.States_WF, 2) != length(wf.Names_WF)
        throw(DimensionMismatch("Number of wind field states in States_WF ($(size(wf.States_WF, 1))) must match length of Names_WF ($(length(wf.Names_WF)))"))
    end
    
    # Create DataFrame with wind field state names as columns
    df = DataFrame()
    
    # Add each wind field state as a column
    for (i, wf_name) in enumerate(wf.Names_WF)
        df[!, Symbol(wf_name)] = wf.States_WF[:, i]
    end
    
    return df
end

"""
    Base.getproperty(wf::WindFarm, name::Symbol)

Custom property accessor for WindFarm that provides special properties like `turbines` and `windfield`.

# Special Properties
- `wf.turbines`: Returns a DataFrame with turbine names as columns and turbine states as data
- `wf.windfield`: Returns a DataFrame with wind field state names as columns and wind field data as data

# Example
```julia
julia> wf.turbines
1800×3 DataFrame
  Row │ a        yaw      TI       
      │ Float64  Float64  Float64  
──────┼────────────────────────────
    1 │    0.33      0.0  0.101532
    2 │    0.33      0.0  0.101532
    3 │    0.33      0.0  0.101532
  ⋮  │    ⋮        ⋮        ⋮ 

julia> wf.windfield
1800×4 DataFrame
  Row │ wind_vel  wind_dir  TI0      OP_ori  
      │ Float64   Float64   Float64  Float64 
──────┼──────────────────────────────────────
    1 │      8.2     195.0    0.062    195.0
    2 │      8.2     195.0    0.062    195.0
  ⋮  │      ⋮          ⋮        ⋮       ⋮ 
```
"""
function Base.getproperty(wf::WindFarm, name::Symbol)
    if name === :turbines
        return turbines(wf)
    elseif name === :windfield
        return windfield(wf)
    else
        # Fall back to default behavior for all other properties
        return getfield(wf, name)
    end
end
