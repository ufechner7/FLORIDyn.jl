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
- `phi`: Wind direction values (in degrees) for all turbines at the specified time

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
        # Use a range (1:wf.nT) to avoid allocating an index vector
        phi = getWindDirT(set.dir_mode, wind.dir, 1:wf.nT, t)
    end

    return phi
end

"""
    correctDir!(::Direction_All, set, wf, wind, t)

Apply direction correction to all operational points in the wind farm using the Direction_All strategy.

# Arguments
- `::Direction_All`: The direction correction strategy type that applies corrections to all operational points
- `set::Settings`: Simulation settings containing the direction mode configuration
- `wf::WindFarm`: Wind farm object containing turbine states and configuration (modified in-place)
- `wind::Wind`: Wind field data structure containing direction information and input type
- `t`: Current simulation time for temporal interpolation

# Returns
- `nothing`: This function modifies the wind farm state in-place and returns nothing

# Description
This function applies a direction correction using the Direction_All strategy, which updates the 
wind direction for all operational points in the wind farm. The function performs the following operations:

1. **Data Retrieval**: Calls `getDataDir` to obtain current wind direction data for all turbines
2. **Universal State Update**: Updates the wind direction for all operational points (`wf.States_WF[:, 2]`)
3. **Observation Point Orientation**: If the state matrix has 4 columns, updates the observation 
   point orientation (`wf.States_WF[wf.StartI, 4]`) for turbine starting positions only

The correction applies a uniform direction (`phi[1]`) to all operational points, representing 
a scenario where the entire wind field experiences the same directional conditions.

# Key Behavior
- **All OPs Updated**: Every operational point receives the same direction value
- **Uniform Direction**: Uses `phi[1]` for all operational points regardless of turbine
- **Selective OP Orientation**: Only turbine starting positions get OP orientation updates
- **Complete Coverage**: Affects the entire wind farm state matrix at once

This strategy is suitable for:
- Scenarios with uniform wind direction across the entire wind farm
- Rapid directional changes affecting all operational points simultaneously
- Simplified modeling where spatial direction variations are negligible

# Examples
```julia
# Apply direction correction to all operational points
correctDir!(Direction_All(), settings, wind_farm, wind_data, 100.0)

# All operational points now have identical wind direction
@assert all(wind_farm.States_WF[:, 2] .== phi[1])

# OP orientations match wind direction for turbine starting positions only
if size(wind_farm.States_WF, 2) == 4
    for i in 1:wind_farm.nT
        start_idx = wind_farm.StartI[i, 1]
        @assert wind_farm.States_WF[start_idx, 4] == wind_farm.States_WF[start_idx, 2]
    end
end
```

# Implementation Notes
- Uses first element of direction vector (`phi[1]`) for all operational points
- Broadcasts the same direction value to all rows of column 2
- OP orientation (column 4) is updated only for turbine starting positions (`wf.StartI`)
- More comprehensive than `Direction_None` which only affects starting positions
- Compatible with different wind input modes (constant, interpolated, random walk)

# See also
- [`getDataDir`](@ref): Function for retrieving wind direction data
- [`Direction_None`](@ref), [`Direction_Influence`](@ref): Alternative correction strategies
- [`Settings`](@ref): Simulation settings structure
- [`WindFarm`](@ref): Wind farm configuration structure
- [`Wind`](@ref): Wind field data structure
"""
function correctDir!(::Direction_All, set, wf, wind, t)
    # Get Data
    phi = getDataDir(set, wind, wf, t)
    # Correct
    wf.States_WF[:, 2] .= phi[1]
    # OP Orientation = turbine wind direction
    if size(wf.States_WF, 2) == 4
       wf.States_WF[wf.StartI, 4] .= wf.States_WF[wf.StartI, 2]
    end
    return nothing
end

"""
    correctDir!(::Direction_None, set, wf, wind, t)

Apply wind direction correction only to turbine starting positions without wake-based corrections.

# Arguments
- `::Direction_None`: Direction correction strategy that applies minimal corrections
- `set::Settings`: Simulation settings containing direction mode and interpolation parameters
- `wf::WindFarm`: Wind farm object with turbine states and configuration (modified in-place)
- `wind::Wind`: Wind field data containing direction time series or model parameters
- `t`: Current simulation time for temporal interpolation

# Returns
- `nothing`: Modifies wind farm state in-place

# Description
This correction strategy applies wind direction updates only to the starting positions of each 
turbine (indexed by `wf.StartI`) rather than to all operational points. This represents a 
minimal correction approach that updates turbine base states without affecting downstream 
operational points.

The function performs these operations:
1. **Data Retrieval**: Obtains current wind direction via [`getDataDir`](@ref)
2. **Selective Assignment**: Sets only turbine starting positions (`wf.StartI`) to `phi[1]`
3. **Observation Point Sync**: If present, aligns OP orientation with the uniform direction value

This approach differs from [`Direction_All`](@ref) by only updating turbine starting positions 
rather than all operational points, making it suitable for:
- Scenarios where only turbine base states need direction updates
- Reduced computational overhead compared to full state corrections
- Initial turbine positioning without affecting wake propagation states

# Examples
```julia
# Apply direction correction to turbine starting positions only
correctDir!(Direction_None(), settings, wind_farm, wind_data, 100.0)

# Only starting positions are updated with uniform direction
for i in 1:wind_farm.nT
    start_idx = wind_farm.StartI[i, 1]
    @assert wind_farm.States_WF[start_idx, 2] == phi[1]
end

# OP orientations match the uniform direction (if 4-column state matrix)
if size(wf.States_WF, 2) == 4
    @assert all(wf.States_WF[wf.StartI, 4] .== phi[1])
end
```

# Implementation Notes
- Uses first element of direction vector (`phi[1]`) for all turbine starting positions
- Only affects turbine starting positions (`wf.StartI`), not all operational points
- OP orientation is set to the raw direction value (`phi[1]`), not the turbine state
- Compatible with different wind input modes (constant, interpolated, random walk)
- More computationally efficient than full state corrections

# See also
- [`correctDir!`](@ref): Generic direction correction interface
- [`getDataDir`](@ref): Wind direction data retrieval
- [`Direction_All`](@ref), [`Direction_Influence`](@ref): Alternative correction strategies
"""
function correctDir!(::Direction_None, set, wf, wind, t)
    # Get Data
    phi = getDataDir(set, wind, wf, t)
    # Correct
    wf.States_WF[wf.StartI, 2] .= phi[1]
    # OP Orientation = turbine wind direction
    if size(wf.States_WF, 2) == 4
       wf.States_WF[wf.StartI, 4] .= phi[1]
    end
    return nothing
end

"""
    correctDir!(::Direction_Influence, set, wf, wind, t)

Influence-based wind direction correction. Each turbine's direction may depend on
upstream operating point (OP) combinations defined in `wf.intOPs` and weighted by
`wf.Weight` according to dependency lists in `wf.dep`.

# Behavior per turbine iT
1. No dependencies (`wf.dep[iT]` empty): assign raw interpolated direction `phi[iT]`.
2. Single interpolation row (`size(wf.intOPs[iT],1)==1` with 4 columns): treat row as
   `[idx1 w1 idx2 w2]` and set direction to `w1*dir[idx1] + w2*dir[idx2]` using current
   (already corrected) turbine directions.
3. Multiple rows (each 4 columns) with non-zero total weight: first compute each row's
   local interpolated direction, then take a weight-averaged sum using `wf.Weight[iT]`.
   Zero or missing weights => fallback to raw `phi[iT]`.
4. Malformed / missing interpolation data: fallback to raw `phi[iT]`.

Column 4 (OP orientation) is synchronized to the updated turbine direction when present.

# Arguments
- `::Direction_Influence`: Strategy marker
- `set::Settings`: Simulation settings (used to fetch direction data)
- `wf::WindFarm`: Wind farm state (modified in-place)
- `wind::Wind`: Wind data source
- `t`: Simulation time

# Returns
- `nothing` (updates `wf` in-place)

# Notes
- Assumes rows of `wf.intOPs[iT]` are of the form `[idx1 w1 idx2 w2]`.
- Uses already updated `wf.States_WF[idx,2]` values within the current loop iteration.
- Robust to missing `dep`, `Weight`, or `intOPs` entries (falls back gracefully).
"""
function correctDir!(::Direction_Influence, set::Settings, wf::WindFarm, wind::Wind, t)
    # Base directions for all turbines
    phi = getDataDir(set, wind, wf, t)

    nT = wf.nT
    has_dep = !isempty(wf.dep)
    has_intOPs = !isempty(wf.intOPs)
    has_weights = !isempty(wf.Weight)
    has_orientation = size(wf.States_WF, 2) == 4

    for iT in 1:nT
        dep_i = (has_dep && length(wf.dep) >= iT) ? wf.dep[iT] : Int[]
        start_idx = wf.StartI[iT, 1]

        if isempty(dep_i)
            # No dependencies -> raw direction
            wf.States_WF[start_idx, 2] = phi[iT]
            if has_orientation
                wf.States_WF[start_idx, 4] = wf.States_WF[start_idx, 2]
            end
            continue
        end

        intOPs_i = (has_intOPs && length(wf.intOPs) >= iT) ? wf.intOPs[iT] : Array{Float64,2}(undef,0,0)
        weights_i = (has_weights && length(wf.Weight) >= iT) ? wf.Weight[iT] : Float64[]

        if size(intOPs_i, 1) == 1 && size(intOPs_i, 2) == 4
            r = intOPs_i
            idx1 = Int(r[1]); w1 = r[2]
            idx2 = Int(r[3]); w2 = r[4]
            wf.States_WF[start_idx, 2] = wf.States_WF[idx1, 2] * w1 + wf.States_WF[idx2, 2] * w2
        elseif size(intOPs_i, 1) > 0 && size(intOPs_i, 2) == 4
            if isempty(weights_i) || sum(weights_i) == 0.0
                wf.States_WF[start_idx, 2] = phi[iT]
            else
                sum_wP = 0.0
                sum_w = 0.0
                nrows = size(intOPs_i, 1)
                for iiT in 1:nrows
                    w = (iiT <= length(weights_i)) ? weights_i[iiT] : 0.0
                    w == 0.0 && continue
                    row = intOPs_i[iiT, :]
                    idx1 = Int(row[1]); w1 = row[2]
                    idx2 = Int(row[3]); w2 = row[4]
                    # Validate indices before accessing States_WF
                    if idx1 > 0 && idx1 <= size(wf.States_WF, 1) && idx2 > 0 && idx2 <= size(wf.States_WF, 1)
                        local_dir = wf.States_WF[idx1, 2] * w1 + wf.States_WF[idx2, 2] * w2
                    else
                        continue  # Skip invalid indices
                    end
                    sum_wP += w * local_dir
                    sum_w += w
                end
                wf.States_WF[start_idx, 2] = sum_w > 0 ? (sum_wP / sum_w) : phi[iT]
            end
        else
            # Fallback for malformed/missing intOPs
            wf.States_WF[start_idx, 2] = phi[iT]
        end

        if has_orientation
            wf.States_WF[start_idx, 4] = wf.States_WF[start_idx, 2]
        end
    end
    return nothing
end

