# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    TurbineArray

A structure representing the configuration and properties of a wind turbine array.

# Fields
- `pos::Matrix{Float64}`: A matrix containing the positions of turbines. Each row represents
  a turbine with columns for x, y, and z coordinates (in meters).
- `type::Vector{String}`: A vector of strings specifying the type/model of each turbine.
- `init_States::Matrix{Float64}`: A matrix containing the initial states of each turbine.
  Each row represents a turbine with columns for:
  - Column 1: `a` - axial induction factor
  - Column 2: `yaw` - initial yaw angle (in degrees)
  - Column 3: `ti` - turbulence intensity

# Example
```julia
# Create a simple 2-turbine array
pos = [0.0 0.0 0.0; 500.0 0.0 0.0]  # Two turbines 500m apart
type = ["NREL_5MW", "NREL_5MW"]
init_states = [0.33 0.0 0.1; 0.33 0.0 0.1]  # Both start with same initial conditions
turbines = TurbineArray(pos, type, init_states)
```

# See Also
- [`turbineArrayProperties`](@ref): Function to load turbine array data from YAML files
"""
struct TurbineArray
    pos::Matrix{Float64}
    type::Vector{String}
    init_States::Matrix{Float64}
end

"""
    turbineArrayProperties(filepath::String) -> TurbineArray

Reads turbine array properties from a YAML configuration file and returns a structured representation.

# Arguments
- `filepath::String`: The path to the `.yaml` file containing turbine array data.

# Returns
- `TurbineArray`: A structured object containing the turbine array configuration with fields:
  - `pos`: Matrix of turbine positions (x, y, z coordinates in meters)
  - `type`: Vector of turbine type/model names
  - `init_States`: Matrix of initial turbine states (axial induction factor, yaw angle, turbulence intensity)

# YAML File Format
The input YAML file should have the following structure:
```yaml
turbines:
  - x: 0.0      # x-coordinate (m)
    y: 0.0      # y-coordinate (m) 
    z: 0.0      # z-coordinate (m)
    type: "NREL_5MW"  # turbine type/model
    a: 0.33     # axial induction factor
    yaw: 0.0    # initial yaw angle (degrees)
    ti: 0.1     # turbulence intensity
  - x: 500.0
    y: 0.0
    z: 0.0
    type: "NREL_5MW"
    a: 0.33
    yaw: 0.0
    ti: 0.1
```

# Example
```julia
# Load turbine configuration from file
turbines = turbineArrayProperties("config/wind_farm.yaml")

# Access turbine positions
println("Number of turbines: ", size(turbines.pos, 1))
println("First turbine position: ", turbines.pos[1, :])

# Access turbine types
println("Turbine types: ", turbines.type)
```

# See Also
- [`TurbineArray`](@ref): The structure returned by this function
"""
function turbineArrayProperties(filepath::String)
    # Load YAML data
    data = YAML.load_file(filepath)

    turbines = data["turbines"]

    # Extract Position: x, y, z
    Pos = [Float64[t["x"], t["y"], t["z"]] for t in turbines]
    Pos = reduce(vcat, [p' for p in Pos])  # transpose and concatenate into matrix (9×3)

    # Extract Type
    Type = [String(t["type"]) for t in turbines]

    # Extract Init States: a, yaw, ti
    Init_States = [Float64[t["a"], t["yaw"], t["ti"]] for t in turbines]
    Init_States = reduce(vcat, [s' for s in Init_States])  # transpose and concatenate

    return TurbineArray(Pos, Type, Init_States)
end
