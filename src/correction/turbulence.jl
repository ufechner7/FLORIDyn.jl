# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    correctTI!(::TI_None, set::Settings, wf, wind, t) -> Nothing

Update turbulence intensity values in the wind farm state matrix without correction.

This function implements the "no correction" strategy for turbulence intensity, where 
the wind farm turbulence intensity values are updated with fresh data from the wind 
field model without applying any correction algorithms. It serves as the baseline 
approach for turbulence intensity handling in FLORIDyn simulations.

# Arguments
- `::TI_None`: Dispatch type indicating no turbulence intensity correction algorithm
- `set::Settings`: Settings object containing simulation configuration and turbulence model parameters
    - `set.turb_mode`: Turbulence model configuration specifying the retrieval method
- `wf::WindFarm`: Wind farm object containing the state matrices to be updated
    - `wf.States_WF`: Wind field states matrix where column 3 contains turbulence intensity values
    - `wf.StartI`: Starting indices for each turbine's observation points
    - `wf.nT`: Number of turbines
- `wind::Wind`: Wind configuration object containing turbulence intensity data
    - `wind.ti`: Turbulence intensity data or model parameters
- `t`: Current simulation time for time-dependent turbulence intensity retrieval

# Returns
- `Nothing`: The function modifies the wind farm state in-place

# Behavior
1. Retrieves current turbulence intensity values for all turbines using [`getDataTI`](@ref)
2. Updates the wind farm state matrix `wf.States_WF` at rows `wf.StartI` and column 3
3. Transposes the turbulence intensity vector to match the matrix structure
4. Provides error handling for matrix update operations

# Example
```julia
# Update turbulence intensity without correction at t=50s
correctTI!(TI_None(), settings, wf, wind, 50.0)

# The wind farm state matrix is now updated with new TI values
current_ti = wf.States_WF[wf.StartI, 3]
```

# Notes
- The function modifies the wind farm object in-place (indicated by the `!` suffix)
- This "no correction" approach provides baseline turbulence intensity without 
  applying wake-induced corrections or measurement-based adjustments
- Error handling ensures graceful failure if matrix dimensions are incompatible

# See Also
- [`getDataTI`](@ref): Function used to retrieve turbulence intensity data
- [`TI_None`](@ref): Dispatch type for no correction strategy
"""
function correctTI!(::TI_None, set::Settings, wf, wind, t)
    # correctTI! updates the turbulent intensity (TI) value inT.States_WF
    # at the rowT.StartI and column 3, using data from getDataTI.

    # Get new turbulent intensity value
    TI = getDataTI(set, wind, wf, t)

    # Update TI at the turbine (first OP) row for each turbine to avoid shape mismatch.
    # Replicating across additional OP columns (if any) is unnecessary for baseline mode.
    if ndims(wf.StartI) == 2 && size(wf.StartI,2) > 1 && size(wf.StartI,1) >= wf.nT
        for iT in 1:wf.nT
            wf.States_WF[wf.StartI[iT,1], 3] = TI[iT]
        end
    else
        # Fallback: treat StartI as a flat collection of indices (row vector or length >= nT)
        # Support both matrix with one row and plain vector
        flat = vec(wf.StartI)
        for iT in 1:min(wf.nT, length(flat))
            wf.States_WF[flat[iT], 3] = TI[iT]
        end
    end

    return nothing
end

"""
    getDataTI(set::Settings, wind::Wind, wf::WindFarm, t) -> Vector

Retrieve turbulence intensity data for all turbines at the current simulation time.

This function obtains turbulence intensity values for all turbines in the wind farm
using the configured turbulence model and wind field data. It serves as a wrapper
around the underlying turbulence intensity retrieval system.

# Arguments
- `set::Settings`: Settings object containing simulation configuration
    - `set.turb_mode`: Turbulence model configuration specifying the retrieval method
- `wind::Wind`: Wind configuration object containing turbulence intensity data
    - `wind.ti`: Turbulence intensity data, parameters, or model configuration
- `wf::WindFarm`: Wind farm object containing turbine information
    - `wf.nT`: Number of turbines in the wind farm
- `t`: Current simulation time for time-dependent turbulence intensity models

# Returns
- `Vector`: Turbulence intensity values for all turbines (dimensionless, typically 0.05-0.25)

# Behavior
The function creates a vector of all turbine indices `[1, 2, ..., nT]` and retrieves
the corresponding turbulence intensity values using the specified turbulence model.
The actual retrieval method depends on the `set.turb_mode` configuration and can include:
- Constant turbulence intensity
- Time-interpolated values from data files
- Turbine-specific interpolation
- Random walk models with covariance

# Example
```julia
# Get turbulence intensity for all turbines at t=100s
TI_values = getDataTI(settings, wind_config, wind_farm, 100.0)
println("TI for turbine 1: ", TI_values[1])
```

# See Also
- [`getWindTiT`](@ref): Underlying function for turbulence intensity retrieval
- [`correctTI!`](@ref): Function that uses this data to update wind farm states
"""
function getDataTI(set::Settings, wind::Wind, wf::WindFarm, t)
    TI = getWindTiT(set.turb_mode, wind.ti, collect(1:wf.nT), t)
    return TI
end

"""
    correctTI!(::TI_Influence, set::Settings, wf::WindFarm, wind::Wind, t)

Apply influence-based turbulence intensity correction where each turbine's turbulence intensity
depends on upstream operational points through dependency relationships and interpolation weights.

# Arguments
- `::TI_Influence`: Turbulence intensity correction strategy that uses dependency-based influence modeling
- `set::Settings`: Simulation settings containing turbulence intensity mode and interpolation parameters
- `wf::WindFarm`: Wind farm object with turbine states and dependency configuration (modified in-place)
- `wind::Wind`: Wind field data containing turbulence intensity time series or model parameters
- `t`: Current simulation time for temporal interpolation

# Returns
- `nothing`: Modifies wind farm state in-place

# Description
This correction strategy implements dependency-based turbulence intensity modeling where each
turbine's ambient turbulence intensity is influenced by upstream operational points. The function
processes each turbine individually based on its dependency configuration and updates column 3
of the wind farm state matrix (`wf.States_WF[:, 3]`).

## Processing Logic per Turbine
1. **No Dependencies**: If `wf.dep[iT]` is empty, assigns the raw ambient turbulence intensity `TI[iT]`
2. **Single Influence Row**: For `wf.intOPs[iT]` with one 4-column row `[idx1 w1 idx2 w2]`,
   computes turbulence intensity as `w1 * States_WF[idx1,3] + w2 * States_WF[idx2,3]`
3. **Multiple Influence Rows**: For multiple 4-column rows, computes the arithmetic mean of
   per-row interpolated turbulence intensities (simple averaging without additional weighting)
4. **Fallback Cases**: Uses raw ambient turbulence intensity `TI[iT]` for malformed data or invalid indices

## Key Features
- **Index Validation**: Validates all operational point indices before accessing state data
- **Robust Handling**: Gracefully handles missing or malformed dependency/interpolation data
- **Sequential Processing**: Uses already-corrected turbine TI values within the same iteration
- **Bounds Checking**: Prevents out-of-bounds access to the wind farm state matrix
- **Simple Averaging**: Multiple influence rows are combined using arithmetic mean (MATLAB-consistent)

## Data Structure Requirements
- `wf.dep[iT]`: Vector of turbine indices that influence turbine `iT`
- `wf.intOPs[iT]`: Matrix where each row is `[idx1 w1 idx2 w2]` for linear interpolation
- `wf.StartI`: 1×nT matrix containing starting row indices for each turbine

# Examples
```julia
# Single influence case - turbine 2 influenced by OPs at indices 1 and 3
wf.dep[2] = [1, 3]
wf.intOPs[2] = [1 0.4 3 0.6]  # 40% from OP 1, 60% from OP 3
correctTI!(TI_Influence(), settings, wf, wind, 100.0)

# Multiple influence case with mean averaging
wf.dep[3] = [1, 2, 4]
wf.intOPs[3] = [1 0.3 2 0.7; 2 0.2 4 0.8]  # Two influence combinations
# Result: mean of (0.3*TI[1] + 0.7*TI[2]) and (0.2*TI[2] + 0.8*TI[4])
correctTI!(TI_Influence(), settings, wf, wind, 100.0)
```

# Implementation Details
- Uses `StartI[1, iT]` indexing (1×nT matrix layout) to access turbine starting positions
- Performs bounds checking on all operational point indices before state access
- Skips invalid indices and counts valid ones for proper averaging in multiple-row cases
- Does not use distance weighting or sophisticated spatial models (matches MATLAB behavior)
- Compatible with different turbulence intensity input modes (constant, interpolated, etc.)

# Error Handling
- Invalid operational point indices are skipped with fallback to ambient turbulence intensity
- Missing dependency/interpolation arrays are treated as empty (no influence)
- Malformed interpolation matrices (wrong dimensions) trigger fallback behavior
- Zero valid interpolations in multiple-row case falls back to ambient TI

# Turbulence Intensity Integration
- Base TI values obtained via `getDataTI` which supports multiple input modes
- Compatible with constant, interpolated, and time-varying turbulence intensity sources
- Turbulence intensity values are dimensionless (typically 0.05-0.25)

# Notes
- Unlike velocity and direction corrections, this function does not use Weight arrays for combining multiple rows
- The MATLAB reference uses simple mean averaging for multiple influence combinations
- Order of turbine processing matters as earlier turbines' corrected TI affects later calculations
- More computationally intensive than `TI_None` due to dependency processing

# See also
- [`getDataTI`](@ref): Function for retrieving ambient turbulence intensity data
- [`correctTI!(::TI_None, set::Settings, wf, wind, t)`](@ref): Simpler TI correction without influence
- [`Settings`](@ref): Simulation settings structure containing turbulence intensity mode
- [`WindFarm`](@ref): Wind farm configuration structure with dependency data
- [`Wind`](@ref): Wind field data structure containing turbulence intensity information
"""
function correctTI!(::TI_Influence, set::Settings, wf::WindFarm, wind::Wind, t)
    # Base TI values for all turbines
    TI = getDataTI(set, wind, wf, t)

    nT = wf.nT
    has_dep = !isempty(wf.dep)
    has_intOPs = !isempty(wf.intOPs)
    
    nStatesWF = size(wf.States_WF, 1)

    for iT in 1:nT
        start_idx = wf.StartI[1, iT]  # FIX: correct indexing (row 1, turbine iT)
        dep_i = (has_dep && length(wf.dep) >= iT) ? wf.dep[iT] : Int[]
        if isempty(dep_i)
            wf.States_WF[start_idx, 3] = TI[iT]
            continue
        end

        intOPs_i = (has_intOPs && length(wf.intOPs) >= iT) ? wf.intOPs[iT] : Array{Float64,2}(undef,0,0)
        nrows = size(intOPs_i, 1)

        if nrows == 1 && size(intOPs_i, 2) == 4
            row = intOPs_i
            idx1 = Int(row[1]); w1 = row[2]
            idx2 = Int(row[3]); w2 = row[4]
            # Validate indices before accessing States_WF
            valid = 0 < idx1 <= nStatesWF && 0 < idx2 <= nStatesWF
            if valid
                wf.States_WF[start_idx, 3] = wf.States_WF[idx1, 3] * w1 + wf.States_WF[idx2, 3] * w2
            else
                wf.States_WF[start_idx, 3] = TI[iT]  # fallback
            end
        elseif nrows > 1 && size(intOPs_i, 2) == 4
            acc = 0.0
            valid_count = 0
            for r in 1:nrows
                row = intOPs_i[r, :]
                idx1 = Int(row[1]); w1 = row[2]
                idx2 = Int(row[3]); w2 = row[4]
                # Validate indices before accessing States_WF
                if 0 < idx1 <= nStatesWF && 0 < idx2 <= nStatesWF
                    acc += wf.States_WF[idx1, 3] * w1 + wf.States_WF[idx2, 3] * w2
                    valid_count += 1
                end
            end
            wf.States_WF[start_idx, 3] = valid_count > 0 ? (acc / valid_count) : TI[iT]
        else
            # Fallback to base TI if malformed or missing
            wf.States_WF[start_idx, 3] = TI[iT]
        end
    end
    return nothing
end

# function T = correctTi(T,Wind,SimTime)
# %CORRECTTI Correction of the turbulent intensity
# %% Get Data
# TI = getDataTI(Wind,T,SimTime);
# %% Correct
# for iT = 1:T.nT
#     if isempty(T.dep{iT})
#         T.States_WF(T.StartI(iT),3) = TI(iT);
#         continue
#     end
    
#     % Assign amb. TI based on what the influencing OPs carry
#     if size(T.intOPs{iT},1)==1
#         % Only one turbine influencing
#         T.States_WF(T.StartI(iT),3) = ...
#             T.States_WF(T.intOPs{iT}(1),3) * T.intOPs{iT}(2) + ...
#             T.States_WF(T.intOPs{iT}(3),3) * T.intOPs{iT}(4);
#     else
#         % More than one, currently only taking the mean, but could be more
#         % sophisticated (e.g. only from uninfluenced wind turbines,
#         % weighted by distance to OP)
#         T.States_WF(T.StartI(iT),1) = mean(...
#             T.States_WF(T.intOPs{iT}(:,1),3) .* T.intOPs{iT}(:,2) + ...
#             T.States_WF(T.intOPs{iT}(:,3),3) .* T.intOPs{iT}(:,4));
#     end
# end
# end

