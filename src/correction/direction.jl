# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getDataDir(set::Settings, wind::Wind, wf::WindFarm, t)

Retrieve wind direction data for all turbines at the current simulation time.

# Arguments
- `set::Settings`: Simulation settings containing the direction mode configuration
- `wind::Wind`: Wind field data structure containing direction information and input type
- `wf::WindFarm`: Wind farm object containing turbine states and configuration
- `t`: Current simulation time for temporal interpolation

# Returns
- `phi`: Wind direction values (typically in radians) for all turbines at the specified time

# Description
This function reads wind direction data and returns the current wind direction angle (phi) for all 
turbines in the wind farm. The function handles two different input modes:

1. **Random Walk with Mean mode** (`wind.input_dir == "RW_with_Mean"`): Uses the current wind farm 
   state from `wf.States_WF[wf.StartI, 2]` along with the wind direction data to compute direction.

2. **Standard temporal interpolation mode**: Uses the wind direction data directly with temporal 
   interpolation for all turbines at the specified simulation time.

The function dispatches to `getWindDirT` with appropriate parameters based on the input mode, 
ensuring consistent wind direction estimation across different modeling approaches.

# Examples
```julia
# Get wind direction for all turbines at current simulation time
phi = getDataDir(settings, wind_data, wind_farm, 100.0)

# The returned phi contains direction values for all turbines
direction_turbine_1 = phi[1]
```

# Notes
- The function automatically handles different wind input modes through conditional logic
- For random walk mode, uses existing wind farm state as reference
- For standard mode, performs temporal interpolation across all turbines

# See also
- [`getWindDirT`](@ref): Underlying function for wind direction temporal interpolation
- [`Settings`](@ref): Configuration structure containing direction mode settings
- [`Wind`](@ref): Wind field data structure containing direction information and input type
- [`WindFarm`](@ref): Wind farm configuration structure
"""
function getDataDir(set::Settings, wind::Wind, wf::WindFarm, t)
    # Reads wind data and returns the current phi for all turbines

    if wind.input_dir == "RW_with_Mean"
        phi = getWindDirT(set.dir_mode,wf.States_WF[wf.StartI, 2], wind.dir)
    else
        phi = getWindDirT(set.dir_mode, wind.dir, collect(1:wf.nT), t)
    end

    return phi
end

"""
    correctDir!(::Direction_All, set::Settings, wf::WindFarm, wind::Wind, t)

Apply direction correction to all turbines in the wind farm using the Direction_All strategy.

# Arguments
- `::Direction_All`: The direction correction strategy type that applies corrections to all turbines
- `set::Settings`: Simulation settings containing the direction mode configuration
- `wf::WindFarm`: Wind farm object containing turbine states and configuration (modified in-place)
- `wind::Wind`: Wind field data structure containing direction information and input type
- `t`: Current simulation time for temporal interpolation

# Returns
- `nothing`: This function modifies the wind farm state in-place and returns nothing

# Description
This function applies a direction correction using the Direction_All strategy, which updates the 
wind direction for all turbines in the wind farm. The function performs the following operations:

1. **Data Retrieval**: Calls `getDataDir` to obtain current wind direction data for all turbines
2. **State Update**: Updates the wind direction in the wind farm state (`wf.States_WF[:, 2]`)
3. **Operational Point Orientation**: If the state matrix has 4 columns, also updates the 
   operational point orientation (`wf.States_WF[wf.StartI, 4]`) to match the wind direction

The correction is applied uniformly to all turbines using the first direction value from the 
retrieved direction data.

# Examples
```julia
# Apply direction correction to all turbines
correctDir!(Direction_All(), settings, wind_farm, wind_data, 100.0)

# The wind farm state is modified in-place
current_direction = wind_farm.States_WF[1, 2]  # Updated direction for first turbine
```

# Notes
- This function modifies the wind farm state in-place (indicated by the `!` suffix)
- All turbines receive the same direction correction value (`phi[1]`)
- The operational point orientation is only updated if the state matrix has 4 columns
- Direction values are typically in radians following standard wind engineering conventions

# See also
- [`getDataDir`](@ref): Function for retrieving wind direction data
- [`Direction_All`](@ref): Direction correction strategy type
- [`Settings`](@ref): Simulation settings structure
- [`WindFarm`](@ref): Wind farm configuration structure
- [`Wind`](@ref): Wind field data structure
"""
function correctDir!(::Direction_All, set::Settings, wf::WindFarm, wind::Wind, t)
    # Get Data
    phi = getDataDir(set, wind, wf, t)
    # Correct
    wf.States_WF[:, 2] .= phi[1]
    # OP Orientation = turbine wind direction
    if size(wf.States_WF, 2) == 4
       wf.States_WF[wf.StartI, 4] .= phi[1]
    end
    return nothing
end

"""
    correctDir!(ub::UnifiedBuffers, ::Direction_All, set::Settings, wf::WindFarm, wind::Wind, t)

Buffered version of correctDir! that uses pre-allocated wind direction buffer.
"""
function correctDir!(ub::UnifiedBuffers, ::Direction_All, set::Settings, wf::WindFarm, wind::Wind, t)
    # Get Data using buffered version
    phi = getDataDir!(ub.wind_dir_buffer, set, wind, wf, t)
    # Correct
    wf.States_WF[:, 2] .= phi[1]
    # OP Orientation = turbine wind direction
    if size(wf.States_WF, 2) == 4
       wf.States_WF[wf.StartI, 4] .= phi[1]
    end
    return nothing
end

"""
    getDataDir!(buffer::Vector{Float64}, set::Settings, wind::Wind, wf::WindFarm, t)

Buffered version of getDataDir that uses pre-allocated buffer for output.
"""
function getDataDir!(buffer::Vector{Float64}, set::Settings, wind::Wind, wf::WindFarm, t)
    # Buffered version that reuses pre-allocated buffer
    if wind.input_dir == "RW_with_Mean"
        # For RW_with_Mean mode, use regular function and copy to buffer
        phi = getWindDirT(set.dir_mode, wf.States_WF[wf.StartI, 2], wind.dir)
        if isa(phi, Number)
            fill!(buffer, phi)
        else
            buffer[1:length(phi)] .= phi
        end
        return buffer
    else
        # For interpolation mode, use buffered version
        return getWindDirT_buffered!(buffer, set.dir_mode, wind.dir, collect(1:wf.nT), t)
    end
end
