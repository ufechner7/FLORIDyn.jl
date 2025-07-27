# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getDataDir(set::Settings, wind::Wind, wf, t)

Retrieve wind direction data for all turbines at the current simulation time.

# Arguments
- `set::Settings`: Simulation settings containing the direction mode configuration
- `wind::Wind`: Wind field data structure containing direction information and input type
- `wf`: Wind farm object containing turbine states
- `t`: Current simulation time for temporal interpolation [s]

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
phi = getDataDir(settings, wind, wf, 100.0)

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
"""
function getDataDir(set::Settings, wind::Wind, wf, t)
    # Reads wind data and returns the current phi for all turbines

    if wind.input_dir == "RW_with_Mean"
        phi = getWindDirT(set.dir_mode,wf.States_WF[wf.StartI, 2], wind.dir)
    else
        phi = getWindDirT(set.dir_mode, wind.dir, collect(1:wf.nT), t)
    end

    return phi
end

"""
    correctDir!(::Direction_All, set::Settings, wf, Wind, t)

Corrects the direction based on the provided parameters.

# Arguments
- `::Direction_All`: The direction correction strategy or type.
- `set`:     The settings for the simulation.`
- `wf`:      The [WindFarm](@ref)
- `wind`:    The wind data or wind state.
- `t`:       The simulation time.

# Description
This function applies a direction correction using the specified strategy, updating the state in-place.
"""
function correctDir!(::Direction_All, set::Settings, wf, wind, t)
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

# function correctDir!(::Direction_All, set::Settings, wf, Wind, SimTime)
#     # CORRECTDIR Correction of the wind direction

#     ## Get Data
#     phi = getDataDir(set, Wind, wf, SimTime)

#     ## Correct the wind direction in the turbine states
#     wf.States_WF[:, 2] .= phi[1]

#     # OP Orientation = turbine wind direction
#     if size(wf.States_WF, 2) == 4
#        wf.States_WF[wf.StartI, 4] .= phi[1]
#     end

#     return wf
# end

