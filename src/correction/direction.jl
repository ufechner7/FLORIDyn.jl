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
        # Use a range (1:wf.nT) to avoid allocating an index vector
        phi = getWindDirT(set.dir_mode, wind.dir, 1:wf.nT, t)
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
3. **Observation Point Orientation**: If the state matrix has 4 columns, also updates the 
   observation point orientation (`wf.States_WF[wf.StartI, 4]`) to match the wind direction

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
- The observation point orientation is only updated if the state matrix has 4 columns
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
    correctDir!(::Direction_None, set::Settings, wf::WindFarm, wind::Wind, t)

Apply basic direction correction using the Direction_None strategy.

# Arguments
- `::Direction_None`: The direction correction strategy type that applies minimal corrections
- `set::Settings`: Simulation settings containing the direction mode configuration
- `wf::WindFarm`: Wind farm object containing turbine states and configuration (modified in-place)
- `wind::Wind`: Wind field data structure containing direction information and input type
- `t`: Current simulation time for temporal interpolation

# Returns
- `nothing`: This function modifies the wind farm state in-place and returns nothing

# Description
This function applies a basic direction correction using the Direction_None strategy. Despite 
its name suggesting "no correction," this function still performs essential direction updates:

1. **Data Retrieval**: Calls `getDataDir` to obtain current wind direction data for all turbines
2. **State Update**: Updates the wind direction for all turbines (`wf.States_WF[:, 2]`) to the 
   first direction value from the retrieved data
3. **Observation Point Orientation**: If the state matrix has 4 columns, sets the observation 
   point orientation (`wf.States_WF[wf.StartI, 4]`) to match the current turbine wind direction

The key difference from `Direction_All` is in the observation point orientation handling: 
instead of using the retrieved direction data directly, it uses the current turbine state.

# Examples
```julia
# Apply basic direction correction
correctDir!(Direction_None(), settings, wind_farm, wind_data, 100.0)

# Check the updated state
turbine_direction = wind_farm.States_WF[1, 2]  # All turbines have same direction
op_orientation = wind_farm.States_WF[wind_farm.StartI, 4]  # Matches turbine direction
```

# Notes
- This function modifies the wind farm state in-place (indicated by the `!` suffix)
- All turbines receive the same direction value (`phi[1]`)
- The observation point orientation copies the turbine direction rather than using raw data
- Direction values are typically in radians following standard wind engineering conventions
- "None" refers to minimal correction strategy, not absence of all direction updates

# See also
- [`correctDir!(::Direction_All, ...)`](@ref): Alternative direction correction strategy
- [`getDataDir`](@ref): Function for retrieving wind direction data
- [`Direction_None`](@ref): Direction correction strategy type
- [`Settings`](@ref): Simulation settings structure
- [`WindFarm`](@ref): Wind farm configuration structure
- [`Wind`](@ref): Wind field data structure
"""
function correctDir!(::Direction_None, set::Settings, wf::WindFarm, wind::Wind, t)
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
    correctDir!(::Direction_Influence, set::Settings, wf::WindFarm, wind::Wind, t)

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

    @inbounds for iT in 1:nT
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
                @inbounds for iiT in 1:nrows
                    w = (iiT <= length(weights_i)) ? weights_i[iiT] : 0.0
                    w == 0.0 && continue
                    row = intOPs_i[iiT, :]
                    idx1 = Int(row[1]); w1 = row[2]
                    idx2 = Int(row[3]); w2 = row[4]
                    local_dir = wf.States_WF[idx1, 2] * w1 + wf.States_WF[idx2, 2] * w2
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

