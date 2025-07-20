# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    turbineArrayProperties(filepath::String)

Reads turbine array properties from the specified `.yaml` file.

# Arguments
- `filepath::String`: The path to the `.yaml` file containing turbine array data.

# Returns
Returns the tuple (Pos, Type, Init_States) of the turbine array as read from the file.

# Notes
Reads a YAML file containing turbine configuration and returns a named tuple
with fields `Pos`, `Type`, and `Init_States`, matching the structure
of the original `turbineArrayProperties()` function.
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

    return (
        Pos = Pos,
        Type = Type,
        Init_States = Init_States
    )
end
