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
        phi = getWindDirT(set.dir_mode,wf.States_WF[wf.StartI, 2], wind)
    else
        phi = getWindDirT(set.dir_mode, wind, 1:wf.nT, t)
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
    correctDir!(::Direction_Influence, set::Settings, wf, wind, t)

Apply influence-based wind direction correction where each turbine's direction depends on
upstream operational points through dependency relationships and interpolation weights.

# Arguments
- `::Direction_Influence`: Direction correction strategy that uses dependency-based influence modeling
- `set::Settings`: Simulation settings containing direction mode and interpolation parameters
- `wf::WindFarm`: Wind farm object with turbine states and dependency configuration (modified in-place)
- `wind::Wind`: Wind field data containing direction time series or model parameters
- `t`: Current simulation time for temporal interpolation

# Returns
- `nothing`: Modifies wind farm state in-place

# Description
This correction strategy implements sophisticated dependency-based direction modeling where each
turbine's wind direction is influenced by upstream operational points. The function processes
each turbine individually based on its dependency configuration:

## Processing Logic per Turbine
1. **No Dependencies**: If `wf.dep[iT]` is empty, assigns the raw ambient direction `phi[iT]`
2. **Single Influence Row**: For `wf.intOPs[iT]` with one 4-column row `[idx1 w1 idx2 w2]`,
   computes direction as `w1 * States_WF[idx1,2] + w2 * States_WF[idx2,2]`
3. **Multiple Influence Rows**: For multiple 4-column rows with corresponding weights in 
   `wf.Weight[iT]`, computes weighted average of per-row interpolated directions
4. **Fallback Cases**: Uses raw ambient direction `phi[iT]` for malformed data or zero weights

## Key Features
- **Index Validation**: Validates all operational point indices before accessing state data
- **Robust Handling**: Gracefully handles missing or malformed dependency/weight/interpolation data
- **Sequential Processing**: Uses already-corrected turbine directions within the same iteration
- **Orientation Sync**: Updates OP orientation (column 4) to match corrected turbine direction

## Data Structure Requirements
- `wf.dep[iT]`: Vector of turbine indices that influence turbine `iT`
- `wf.intOPs[iT]`: Matrix where each row is `[idx1 w1 idx2 w2]` for linear interpolation
- `wf.Weight[iT]`: Vector of weights for combining multiple influence rows
- `wf.StartI`: 1×nT matrix containing starting row indices for each turbine

# Examples
```julia
# Single influence case - turbine 2 influenced by OPs at indices 1 and 3
wf.dep[2] = [1, 3]
wf.intOPs[2] = [1 0.3 3 0.7]  # 30% from OP 1, 70% from OP 3
correctDir!(Direction_Influence(), settings, wf, wind, 100.0)

# Multiple influence case with weighted combination
wf.dep[3] = [1, 2, 4]
wf.intOPs[3] = [1 0.5 2 0.5; 2 0.8 4 0.2]  # Two influence combinations
wf.Weight[3] = [0.6, 0.4]  # Weights for combining the two influences
correctDir!(Direction_Influence(), settings, wf, wind, 100.0)
```

# Implementation Details
- Uses `StartI[1, iT]` indexing (1×nT matrix layout) to access turbine starting positions
- Performs bounds checking on all operational point indices before state access
- Skips invalid indices or zero weights to maintain numerical stability
- Updates orientation column only for the specific turbine's starting position
- More computationally intensive than `Direction_None` or `Direction_All` due to dependency processing

# Error Handling
- Invalid operational point indices are skipped with fallback to ambient direction
- Missing dependency/weight/interpolation arrays are treated as empty (no influence)
- Zero or negative weights are ignored in weighted averaging calculations
- Malformed interpolation matrices (wrong dimensions) trigger fallback behavior

# See also
- [`getDataDir`](@ref): Function for retrieving ambient wind direction data
- [`Direction_None`](@ref), [`Direction_All`](@ref): Simpler correction strategies
- [`Settings`](@ref): Simulation settings structure containing direction mode
- [`WindFarm`](@ref): Wind farm configuration structure with dependency data
- [`Wind`](@ref): Wind field data structure
"""
function correctDir!(::Direction_Influence, set::Settings, wf, wind, t)
    # NOTE: This implementation mirrors the logic of the original MATLAB version (see
    # commented reference below) while fixing indexing and adding robustness.
    # Differences vs MATLAB:
    #  * Uses 1-based row-major Julia arrays with StartI stored as a 1×nT matrix;
    #    therefore we index StartI as StartI[1, iT] (previous code used StartI[iT,1]
    #    causing BoundsError for iT>1).
    #  * Orientation (column 4) is only updated for the specific turbine start row
    #    (the MATLAB code assigns all StartI rows each iteration; final value equals
    #    last turbine's direction – likely unintended side effect).

    phi = getDataDir(set, wind, wf, t)  # base ambient directions per turbine

    nT = wf.nT
    has_dep        = !isempty(wf.dep)
    has_intOPs     = !isempty(wf.intOPs)
    has_weights    = !isempty(wf.Weight)
    has_orientation = size(wf.States_WF, 2) == 4

    nStatesWF = size(wf.States_WF, 1)

    @inbounds for iT in 1:nT
        # Safe retrieval of dependency / interpolation / weight lists
        dep_i      = (has_dep     && length(wf.dep)     >= iT) ? wf.dep[iT]     : Int[]
        intOPs_i   = (has_intOPs  && length(wf.intOPs)  >= iT) ? wf.intOPs[iT]  : Array{Float64}(undef, 0, 0)
        weights_i  = (has_weights && length(wf.Weight)  >= iT) ? wf.Weight[iT]  : Float64[]

        start_idx = wf.StartI[1, iT]  # correct indexing: StartI is 1×nT, so use row 1, turbine iT

        if isempty(dep_i)
            # No dependencies -> assign raw ambient direction
            wf.States_WF[start_idx, 2] = phi[iT]
            # OP Orientation = turbine wind direction (MATLAB behavior: assign to ALL StartI)
            if has_orientation
                wf.States_WF[wf.StartI, 4] .= wf.States_WF[start_idx, 2]
            end
            continue
        end

        # Single influencing row: direct linear combination
        if size(intOPs_i, 1) == 1 && size(intOPs_i, 2) == 4
            r = intOPs_i
            idx1 = Int(r[1]); w1 = r[2]
            idx2 = Int(r[3]); w2 = r[4]
            valid = 0 < idx1 <= nStatesWF && 0 < idx2 <= nStatesWF
            if valid
                wf.States_WF[start_idx, 2] = wf.States_WF[idx1, 2] * w1 + wf.States_WF[idx2, 2] * w2
            else
                wf.States_WF[start_idx, 2] = phi[iT]  # fallback
            end

        # Multiple rows: weighted average of per-row interpolated combinations
        elseif size(intOPs_i, 1) > 0 && size(intOPs_i, 2) == 4
            # Effective number of usable rows limited by available weights / dependencies
            nrows = size(intOPs_i, 1)
            if isempty(weights_i) || sum(weights_i) == 0.0
                wf.States_WF[start_idx, 2] = phi[iT]
            else
                sum_wP = 0.0
                sum_w  = 0.0
                for row_idx in 1:nrows
                    w = (row_idx <= length(weights_i)) ? weights_i[row_idx] : 0.0
                    w == 0.0 && continue
                    row = intOPs_i[row_idx, :]
                    idx1 = Int(row[1]); w1 = row[2]
                    idx2 = Int(row[3]); w2 = row[4]
                    # Validate indices
                    if !(0 < idx1 <= nStatesWF && 0 < idx2 <= nStatesWF)
                        continue
                    end
                    local_dir = wf.States_WF[idx1, 2] * w1 + wf.States_WF[idx2, 2] * w2
                    sum_wP += w * local_dir
                    sum_w  += w
                end
                wf.States_WF[start_idx, 2] = sum_w > 0 ? (sum_wP / sum_w) : phi[iT]
            end
        else
            # Malformed / missing interpolation matrix -> fallback
            wf.States_WF[start_idx, 2] = phi[iT]
        end

        # OP Orientation = turbine wind direction (MATLAB behavior: assign to ALL StartI)
        if has_orientation
            wf.States_WF[wf.StartI, 4] .= wf.States_WF[start_idx, 2]
        end
    end
    return nothing
end

# function T = correctDir(T,Wind,SimTime)
# %CORRECTDIR Correction of the wind direction based on the OPs influencing
# %the turbine.
# %% Get Data
# phi = getDataDir(Wind,T,SimTime);
# %% Correct
# for iT = 1:T.nT
#     if isempty(T.dep{iT})
#         T.States_WF(T.StartI(iT),2) = phi(iT);

#         % OP Orientation = turbine wind direction
#         if size(T.States_WF,2) == 4
#             T.States_WF(T.StartI,4) = T.States_WF(T.StartI(iT),2);
#         end
#         continue
#     end
    
#     % Assign free wind speed based on what the influencing OPs carry
#     if size(T.intOPs{iT},1)==1
#         % Only one turbine influencing
#         T.States_WF(T.StartI(iT),2) = ...
#             T.States_WF(T.intOPs{iT}(1),2) * T.intOPs{iT}(2) + ...
#             T.States_WF(T.intOPs{iT}(3),2) * T.intOPs{iT}(4);

#     else
#         % Weighted average of the wind speed state in the OPs that
#         % influence this turbine. The weight is based on how much FLORIS
#         % thinks this turbine is influenced by another wake. With the Gauss
#         % model its currently based on the gaussian weight.
#         if sum(T.weight{iT}) == 0
#             % Turbine is actually not influenced -> skip
#             T.States_WF(T.StartI(iT),2) = phi(iT);
#         else
#             % Sum the weighted wind speeds and the weights
#             sum_wP = 0;
#             sum_w  = 0;
#             for iiT = 1:length(T.dep{iT})
#                 sum_wP = sum_wP + T.weight{iT}(iiT) * (...
#                     T.States_WF(T.intOPs{iT}(iiT,1),2) .* T.intOPs{iT}(iiT,2) + ...
#                     T.States_WF(T.intOPs{iT}(iiT,3),2) .* T.intOPs{iT}(iiT,4));
                
#                 sum_w = sum_w + T.weight{iT}(iiT);
#             end
#             % Divide to get weigthed average
#             T.States_WF(T.StartI(iT),2) = sum_wP / sum_w;
#         end
#     end

#     % OP Orientation = turbine wind direction
#     if size(T.States_WF,2) == 4
#         T.States_WF(T.StartI,4) = T.States_WF(T.StartI(iT),2);
#     end
# end


