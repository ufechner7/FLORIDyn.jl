# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Parameters

"""
    Measurement

Structure representing a single measurement visualization configuration.

This struct defines the configuration for individual measurement visualizations in wind farm simulations,
including the measurement type and whether it should be displayed in separated plots.

# Fields
- `name::String`: The name/identifier of the measurement type (e.g., "measurements_vel_reduction",
  "measurements_added_turbulence", "measurements_eff_wind_speed")
- `separated::Bool`: Whether the measurement should be plotted in separate individual plots (true) 
  or combined in a single plot (false). Default is false for combined plotting.

# Examples
```julia
# Velocity reduction measurement with separated plots
vel_measurement = Measurement("measurements_vel_reduction", true)

# Turbulence measurement with combined plots  
turb_measurement = Measurement("measurements_added_turbulence", false)
```

# See Also
- [`parse_measurements`](@ref): Function to convert YAML measurement configurations to Measurement arrays
- [`plotMeasurements`](@ref): Function that uses these configurations for visualization
"""
@with_kw struct Measurement
    name::String = ""
    separated::Bool = false
end

# Custom constructor for single string argument
Measurement(name::String) = Measurement(name, false)

"""
    parse_measurements(measurements_yaml::Vector) -> Vector{Measurement}

Convert YAML measurement configurations into an array of Measurement structs.

This function parses the measurement section from a YAML visualization configuration and converts
it into a structured array of Measurement objects. It handles both simple string entries and
complex dictionary entries with name and separated flag specifications.

# Arguments
- `measurements_yaml::Vector`: Vector containing measurement configurations from YAML.
  Each element can be either:
  - A simple string (measurement name)
  - A dictionary with "name" and "separated" keys

# Returns
- `Vector{Measurement}`: Array of Measurement structs with parsed configurations

# YAML Format Support
The function supports multiple YAML formats:

## Simple Format (separated=false by default)
```yaml
measurements:
  - "measurements_vel_reduction"
  - "measurements_added_turbulence" 
```

## Extended Format (with separated flag)
```yaml
measurements:
  - name: "measurements_vel_reduction"
    separated: true
  - name: "measurements_added_turbulence"
    separated: false
```

## Mixed Format
```yaml
measurements:
  - "measurements_vel_reduction"  # Simple string, separated=false
  - name: "measurements_added_turbulence"  # Dictionary format
    separated: true
```

# Algorithm
1. **Iterate** through each measurement entry in the YAML vector
2. **Check type**: Determine if entry is string or dictionary
3. **Extract values**: Get name and separated flag based on format
4. **Create struct**: Construct Measurement struct with parsed values

# Error Handling
- Handles missing separated flags by defaulting to false
- Processes mixed simple/complex YAML structures gracefully
- Validates measurement name extraction

# Examples
```julia
# Parse from YAML data
yaml_measurements = [
    "measurements_vel_reduction",
    Dict("name" => "measurements_added_turbulence", "separated" => true)
]

measurements = parse_measurements(yaml_measurements)
# Returns: [Measurement("measurements_vel_reduction", false), 
#          Measurement("measurements_added_turbulence", true)]
```

# Integration
This function integrates with the existing visualization system:
- Used by Vis struct constructor to parse YAML configurations  
- Outputs compatible with plotMeasurements function
- Supports the existing measurement naming conventions

# See Also
- [`Measurement`](@ref): The struct type created by this function
- [`Vis`](@ref): Visualization settings struct that uses these measurements
- [`plotMeasurements`](@ref): Function that consumes these measurement configurations
"""
function parse_measurements(measurements_yaml::Vector)
    measurements = Measurement[]
    
    for item in measurements_yaml
        if isa(item, String)
            # Simple string format - use default separated=false
            push!(measurements, Measurement(item, false))
        elseif isa(item, Dict)
            # Dictionary format with name and separated keys
            name = get(item, "name", "")
            separated = get(item, "separated", false)
            
            if !isempty(name)
                push!(measurements, Measurement(name, separated))
            end
        end
    end
    
    return measurements
end

"""
    parse_measurements(vis_data::Dict) -> Vector{Measurement}

Parse measurements from a complete visualization configuration dictionary.

This is a convenience method that extracts the measurements section from a full
visualization configuration dictionary and converts it to Measurement structs.

# Arguments
- `vis_data::Dict`: Complete visualization configuration dictionary (e.g., from YAML.load_file)

# Returns
- `Vector{Measurement}`: Array of parsed Measurement structs

# Example
```julia
using YAML
vis_data = YAML.load_file("vis_default.yaml")
measurements = parse_measurements(vis_data["vis"])
```
"""
function parse_measurements(vis_data::Dict)
    if haskey(vis_data, "measurements")
        return parse_measurements(vis_data["measurements"])
    else
        return Measurement[]
    end
end
