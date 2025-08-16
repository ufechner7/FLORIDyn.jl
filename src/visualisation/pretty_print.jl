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

"""
    turbines(wf::WindFarm) -> DataFrame

Create a DataFrame from WindFarm turbine state data with operating point and turbine identifiers.

# Arguments
- `wf::WindFarm`: WindFarm object containing turbine state data

# Returns
- `DataFrame`: DataFrame with columns:
  - `:OP`: Operating point number (1 to nOP, repeated for each turbine)
  - `:Turbine`: Turbine number (1 to nT, each repeated nOP times)
  - Additional columns for each turbine state name from `wf.Names_T`

# Data Structure
The DataFrame rows are organized with all operating points for turbine 1, 
followed by all operating points for turbine 2, etc:
- Rows 1 to nOP: Turbine 1, OP 1 to nOP
- Rows nOP+1 to 2×nOP: Turbine 2, OP 1 to nOP  
- ...and so on

# Throws
- `ArgumentError`: If States_T or Names_T are empty
- `DimensionMismatch`: If dimensions don't match

# Example
```julia
julia> df = wf.turbines
100×5 DataFrame
 Row │ OP     Turbine  a        yaw      TI      
     │ Int64  Int64    Float64  Float64  Float64 
─────┼─────────────────────────────────────────
   1 │     1        1     0.33      0.0   0.101
   2 │     2        1     0.33      0.0   0.102
  ⋮  │   ⋮       ⋮       ⋮        ⋮        ⋮
  50 │    50        1     0.33      0.0   0.105
  51 │     1        2     0.33      0.0   0.101
```
"""
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
    
    # Add OP and Turbine identifier columns
    n_rows = size(wf.States_T, 1)
    n_ops = wf.nOP
    n_turbines = wf.nT
    
    # Create OP column: repeats 1:nOP for each turbine
    op_col = repeat(1:n_ops, n_turbines)
    
    # Create Turbine column: each turbine number repeated nOP times
    turbine_col = repeat(1:n_turbines, inner=n_ops)
    
    # Add identifier columns to DataFrame
    df[!, :OP] = op_col[1:n_rows]  # Truncate to actual row count
    df[!, :Turbine] = turbine_col[1:n_rows]  # Truncate to actual row count
    
    # Add each turbine's states as a column
    for (i, turbine_name) in enumerate(wf.Names_T)
        df[!, Symbol(turbine_name)] = wf.States_T[:, i]
    end
    
    return df
end

"""
    turbines(T::Dict) -> DataFrame

Create a DataFrame from a dictionary containing turbine state data with operating point and turbine identifiers.

# Arguments
- `T::Dict`: Dictionary containing turbine data with keys:
  - `"States_T"` or `:States_T`: Matrix of turbine states (nT*nOP × n_state_variables)
  - `"Names_T"` or `:Names_T`: Vector of turbine state names
  - `"nT"` or `:nT`: Number of turbines 
  - `"nOP"` or `:nOP`: Number of operating points

# Returns
- `DataFrame`: DataFrame with columns:
  - `:OP`: Operating point number (1 to nOP, repeated for each turbine)
  - `:Turbine`: Turbine number (1 to nT, each repeated nOP times)
  - Additional columns for each turbine state name from dictionary

# Data Structure
The DataFrame rows are organized with all operating points for turbine 1, 
followed by all operating points for turbine 2, etc:
- Rows 1 to nOP: Turbine 1, OP 1 to nOP
- Rows nOP+1 to 2×nOP: Turbine 2, OP 1 to nOP  
- ...and so on

# Notes
- Both string keys (`"States_T"`) and symbol keys (`:States_T`) are accepted
- The `States_T` matrix should have dimensions (nT*nOP) × n_state_variables
- Input data is already in the correct row order for the long format

# Throws
- `ArgumentError`: If required keys are missing or data is empty
- `DimensionMismatch`: If dimensions don't match expected structure

# Example
```julia
# Example with MAT file data format
T = Dict(
    "States_T" => [0.33 0.0 0.06; 0.33 0.0 0.06; ...],  # (nT*nOP) × 3 matrix
    "Names_T" => ["a", "yaw", "TI"],                      # State variable names
    "nT" => 3,                                            # Number of turbines  
    "nOP" => 2                                            # Number of operating points
)
df = turbines(T)
# Returns a 6×5 DataFrame with columns: OP, Turbine, a, yaw, TI
```
"""
function turbines(T::Dict)
    # Try to find States_T with different possible key formats
    states_t = nothing
    if haskey(T, "States_T")
        states_t = T["States_T"]
    elseif haskey(T, :States_T)
        states_t = T[:States_T]
    else
        throw(ArgumentError("Dictionary must contain key 'States_T' (or :States_T) with turbine state data"))
    end
    
    # Try to find Names_T with different possible key formats
    names_t = nothing
    if haskey(T, "Names_T")
        names_t = T["Names_T"]
    elseif haskey(T, :Names_T)
        names_t = T[:Names_T]
    else
        throw(ArgumentError("Dictionary must contain key 'Names_T' (or :Names_T) with turbine names"))
    end
    
    # Try to find nT (number of turbines)
    n_turbines = nothing
    if haskey(T, "nT")
        n_turbines = Int(T["nT"])
    elseif haskey(T, :nT)
        n_turbines = Int(T[:nT])
    else
        throw(ArgumentError("Dictionary must contain key 'nT' (or :nT) with number of turbines"))
    end
    
    # Try to find nOP (number of operating points)
    n_ops = nothing
    if haskey(T, "nOP")
        n_ops = Int(T["nOP"])
    elseif haskey(T, :nOP)
        n_ops = Int(T[:nOP])
    else
        throw(ArgumentError("Dictionary must contain key 'nOP' (or :nOP) with number of operating points"))
    end
    
    # Validate inputs
    if isempty(states_t) || isempty(names_t)
        throw(ArgumentError("Dictionary must have non-empty States_T and Names_T to create DataFrame"))
    end
    
    # Check dimensions match
    if size(states_t, 2) != length(names_t)
        throw(DimensionMismatch("Number of state variables in States_T ($(size(states_t, 2))) must match length of Names_T ($(length(names_t)))"))
    end
    
    # Check that States_T has the expected number of rows
    expected_rows = n_turbines * n_ops
    if size(states_t, 1) != expected_rows
        throw(DimensionMismatch("States_T should have $(expected_rows) rows (nT×nOP = $(n_turbines)×$(n_ops)), but has $(size(states_t, 1)) rows"))
    end
    
    # Create DataFrame
    df = DataFrame()
    
    # Create OP column: repeats 1:nOP for each turbine
    op_col = repeat(1:n_ops, n_turbines)
    
    # Create Turbine column: each turbine number repeated nOP times
    turbine_col = repeat(1:n_turbines, inner=n_ops)
    
    # Add identifier columns to DataFrame
    df[!, :OP] = op_col
    df[!, :Turbine] = turbine_col
    
    # Add each state variable as a column (data is already in the correct order)
    for (i, state_name) in enumerate(names_t)
        df[!, Symbol(state_name)] = states_t[:, i]
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

function ops(wf::WindFarm)
    # Check if we have operating points data
    if isempty(wf.States_OP) || isempty(wf.Names_OP)
        throw(ArgumentError("WindFarm must have non-empty States_OP and Names_OP to create DataFrame"))
    end
    
    # Check dimensions match - States_OP should be a column vector or matrix with as many rows as operating points
    # and Names_OP should match the number of columns if States_OP has multiple columns
    if ndims(wf.States_OP) == 2 && size(wf.States_OP, 2) != length(wf.Names_OP)
        throw(DimensionMismatch("Number of operating point variables in States_OP ($(size(wf.States_OP, 2))) must match length of Names_OP ($(length(wf.Names_OP)))"))
    elseif ndims(wf.States_OP) == 1 && length(wf.Names_OP) != 1
        throw(DimensionMismatch("States_OP is a vector but Names_OP has $(length(wf.Names_OP)) elements, expected 1"))
    end
    
    # Create DataFrame with operating point names as columns
    df = DataFrame()
    
    # Handle both vector and matrix cases
    if ndims(wf.States_OP) == 1 || size(wf.States_OP, 2) == 1
        # Single column case
        df[!, Symbol(wf.Names_OP[1])] = vec(wf.States_OP)
    else
        # Multiple columns case
        for (i, op_name) in enumerate(wf.Names_OP)
            df[!, Symbol(op_name)] = wf.States_OP[:, i]
        end
    end
    
    return df
end

"""
    Base.getproperty(wf::WindFarm, name::Symbol)

Custom property accessor for WindFarm that provides special properties like `turbines`, `windfield`, and `ops`.

# Special Properties
- `wf.turbines`: Returns a DataFrame with turbine names as columns and turbine states as data
- `wf.windfield`: Returns a DataFrame with wind field state names as columns and wind field data as data
- `wf.ops`: Returns a DataFrame with operating point names as columns and operating point data as data

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

julia> wf.ops
8760×3 DataFrame
  Row │ U        Dir      TI      
      │ Float64  Float64  Float64 
──────┼──────────────────────────
    1 │     8.5     180.0    0.08
    2 │     9.2     185.5    0.09
  ⋮  │     ⋮         ⋮        ⋮
```
"""
function Base.getproperty(wf::WindFarm, name::Symbol)
    if name === :turbines
        return turbines(wf)
    elseif name === :windfield
        return windfield(wf)
    elseif name === :ops
        return ops(wf)
    else
        # Fall back to default behavior for all other properties
        return getfield(wf, name)
    end
end
