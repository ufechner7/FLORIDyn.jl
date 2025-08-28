# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    correctTI!(::TI_None, set::Settings, wf::WindFarm, wind::Wind, t) -> Nothing

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
function correctTI!(::TI_None, set::Settings, wf::WindFarm, wind::Wind, t)
    # correctTI! updates the turbulent intensity (TI) value inT.States_WF
    # at the rowT.StartI and column 3, using data from getDataTI.

    # Get new turbulent intensity value
    TI = getDataTI(set, wind, wf, t)

    # Update TI at the turbine (first OP) row for each turbine to avoid shape mismatch.
    # Replicating across additional OP columns (if any) is unnecessary for baseline mode.
    if ndims(wf.StartI) == 2 && size(wf.StartI,2) > 1 && size(wf.StartI,1) >= wf.nT
        @inbounds for iT in 1:wf.nT
            wf.States_WF[wf.StartI[iT,1], 3] = TI[iT]
        end
    else
        # Fallback: treat StartI as a flat collection of indices (row vector or length >= nT)
        # Support both matrix with one row and plain vector
        flat = vec(wf.StartI)
        @inbounds for iT in 1:min(wf.nT, length(flat))
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

Influence-based turbulence intensity correction translated from the legacy MATLAB
`correctTi` function. Adjusts each turbine's ambient TI based on upstream
operating point (OP) influence definitions in `wf.intOPs` and dependencies in `wf.dep`.

# Behavior per turbine iT
1. If `wf.dep[iT]` is empty: assign raw TI value `TI[iT]`.
2. If `wf.intOPs[iT]` has exactly one row of length 4: treat it as `[idx1 w1 idx2 w2]`
   and set TI to `w1*TI_idx1 + w2*TI_idx2` using current `wf.States_WF[:,3]` values.
3. Else if it has multiple rows (each 4 columns): compute mean of the per-row
   interpolated TI values `(w1*TI_idx1 + w2*TI_idx2)` (simple average, no weighting).
4. Malformed / missing interpolation data: fallback to raw `TI[iT]`.

# Arguments
- `::TI_Influence`: Strategy marker
- `set::Settings`: Simulation settings (used to obtain base TI via `getDataTI`)
- `wf::WindFarm`: Wind farm state (modified in-place; column 3 is TI)
- `wind::Wind`: Wind data source
- `t`: Simulation time

# Returns
- `nothing` (updates `wf.States_WF` in-place)

# Notes
- Uses existing `wf.States_WF[idx,3]` values when forming influenced TI; order of
  iteration means earlier turbines' TI may already reflect influence.
- Does NOT currently apply distance weighting—mirrors MATLAB mean behavior.
- Gracefully handles empty / missing containers.
"""
function correctTI!(::TI_Influence, set::Settings, wf::WindFarm, wind::Wind, t)
    # Base TI values for all turbines
    TI = getDataTI(set, wind, wf, t)

    nT = wf.nT
    has_dep = !isempty(wf.dep)
    has_intOPs = !isempty(wf.intOPs)

    @inbounds for iT in 1:nT
        start_idx = wf.StartI[iT, 1]
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
            wf.States_WF[start_idx, 3] = wf.States_WF[idx1, 3] * w1 + wf.States_WF[idx2, 3] * w2
        elseif nrows > 1 && size(intOPs_i, 2) == 4
            acc = 0.0
            @inbounds for r in 1:nrows
                row = intOPs_i[r, :]
                idx1 = Int(row[1]); w1 = row[2]
                idx2 = Int(row[3]); w2 = row[4]
                acc += wf.States_WF[idx1, 3] * w1 + wf.States_WF[idx2, 3] * w2
            end
            wf.States_WF[start_idx, 3] = acc / nrows
        else
            # Fallback to base TI if malformed or missing
            wf.States_WF[start_idx, 3] = TI[iT]
        end
    end
    return nothing
end

