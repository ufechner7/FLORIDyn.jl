# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    correctVel(::Velocity_None, set, wf, wind, t, floris, tmpM)

No-op velocity correction: simply overwrites the free-stream wind speed state (first column
of `wf.States_WF`) with the base velocities returned by `getDataVel` for each turbine.

Returns the mutated `wf` and (potentially updated) `wind` for API symmetry with influence
correction. Used when `set.cor_vel_mode isa Velocity_None`.
"""
function correctVel(::Velocity_None, set, wf, Wind, SimTime, paramFLORIS, tmpM)
    # Get data
    U, Wind = getDataVel(set, Wind, wf, SimTime, tmpM, paramFLORIS)

    # Correct Velocity
    wf.States_WF[wf.StartI, 1] = U

    return wf, Wind
end

"""
    correctVel(::Velocity_Influence, set::Settings, wf::WindFarm, wind::Wind, t, floris, tmpM)

Apply influence-based wind velocity correction where each turbine's wind speed depends on
upstream operational points through dependency relationships and interpolation weights.

# Arguments
- `::Velocity_Influence`: Velocity correction strategy that uses dependency-based influence modeling
- `set::Settings`: Simulation settings containing velocity mode and interpolation parameters
- `wf::WindFarm`: Wind farm object with turbine states and dependency configuration (modified in-place)
- `wind::Wind`: Wind field data containing velocity time series or model parameters (may be updated)
- `t`: Current simulation time for temporal interpolation
- `floris`: FLORIS parameter container containing wake model parameters
- `tmpM`: Temporary matrix used in special velocity retrieval modes (e.g., I_and_I estimator)

# Returns
- `(wf, wind)`: Updated wind farm and wind objects (mutates `wf` in-place; returns `wind` for API consistency)

# Description
This correction strategy implements sophisticated dependency-based velocity modeling where each
turbine's ambient wind speed is influenced by upstream operational points. The function processes
each turbine individually based on its dependency configuration and updates column 1 of the
wind farm state matrix (`wf.States_WF[:, 1]`).

## Processing Logic per Turbine
1. **No Dependencies**: If `wf.dep[iT]` is empty, assigns the raw ambient velocity `u[iT]`
2. **Single Influence Row**: For `wf.intOPs[iT]` with one 4-column row `[idx1 w1 idx2 w2]`,
   computes velocity as `w1 * States_WF[idx1,1] + w2 * States_WF[idx2,1]`
3. **Multiple Influence Rows**: For multiple 4-column rows with corresponding weights in 
   `wf.Weight[iT]`, computes weighted average of per-row interpolated velocities
4. **Fallback Cases**: Uses raw ambient velocity `u[iT]` for malformed data or zero weights

## Key Features
- **Index Validation**: Validates all operational point indices before accessing state data
- **Robust Handling**: Gracefully handles missing or malformed dependency/weight/interpolation data

## Data Structure Requirements
- `wf.dep[iT]`: Vector of turbine indices that influence turbine `iT`
- `wf.intOPs[iT]`: Matrix where each row is `[idx1 w1 idx2 w2]` for linear interpolation
- `wf.Weight[iT]`: Vector of weights for combining multiple influence rows
- `wf.StartI`: 1×nT matrix containing starting row indices for each turbine

# Examples
```julia
# Single influence case - turbine 1 influenced by OPs at indices 2 and 4
wf.dep[1] = [2, 4]
wf.intOPs[1] = [2 0.4 4 0.6]  # 40% from OP 2, 60% from OP 4
wf_updated, wind_updated = correctVel(Velocity_Influence(), settings, wf, wind, 100.0, floris, tmpM)

# Multiple influence case with weighted combination
wf.dep[2] = [1, 3, 5]
wf.intOPs[2] = [1 0.3 3 0.7; 3 0.2 5 0.8]  # Two influence combinations
wf.Weight[2] = [0.7, 0.3]  # Weights for combining the two influences
wf_updated, wind_updated = correctVel(Velocity_Influence(), settings, wf, wind, 100.0, floris, tmpM)
```

# Implementation Details
- Uses `StartI[1, iT]` indexing (1×nT matrix layout) to access turbine starting positions
- Performs bounds checking on all operational point indices before state access
- Skips invalid indices or zero weights to maintain numerical stability
- More computationally intensive than `Velocity_None` due to dependency processing
- Calls `getDataVel` to obtain base ambient velocities for all turbines

# Error Handling
- Invalid operational point indices are skipped with fallback to ambient velocity
- Missing dependency/weight/interpolation arrays are treated as empty (no influence)
- Zero or negative weights are ignored in weighted averaging calculations
- Malformed interpolation matrices (wrong dimensions) trigger fallback behavior

# Wind Data Integration
- Base velocities obtained via `getDataVel` which supports multiple input modes
- Compatible with constant, interpolated, EnKF, and estimator-based velocity sources
- May update `wind` object state for certain estimator modes (I_and_I, RW_with_Mean)

# See also
- [`getDataVel`](@ref): Function for retrieving ambient wind velocity data
- [`correctVel(::Velocity_None, ...)`](@ref): Simpler velocity correction without influence
- [`Settings`](@ref): Simulation settings structure containing velocity mode
- [`WindFarm`](@ref): Wind farm configuration structure with dependency data
- [`Wind`](@ref): Wind field data structure
- [`Floris`](@ref): FLORIS parameter structure for wake modeling

# WARNING:
This correction method is not properly tested. Use at your own risk!
"""
function correctVel(::Velocity_Influence, set::Settings, wf::WindFarm, wind::Wind, t, floris, tmpM)
    # Base free wind speeds (may update wind state depending on mode)
    u, wind = getDataVel(set, wind, wf, t, tmpM, floris)

    nT = wf.nT
    has_dep = !isempty(wf.dep)
    has_intOPs = !isempty(wf.intOPs)
    has_weights = !isempty(wf.Weight)
    
    nStatesWF = size(wf.States_WF, 1)

    for iT in 1:nT
        dep_i = (has_dep && length(wf.dep) >= iT) ? wf.dep[iT] : Int[]
        start_idx = wf.StartI[1, iT]  # Indexing: row 1, turbine iT (StartI is a 1×nT matrix)

        if isempty(dep_i)
            wf.States_WF[start_idx, 1] = u[iT]
            continue
        end

        intOPs_i = (has_intOPs && length(wf.intOPs) >= iT) ? wf.intOPs[iT] : Array{Float64,2}(undef,0,0)
        weights_i = (has_weights && length(wf.Weight) >= iT) ? wf.Weight[iT] : Float64[]

        if size(intOPs_i,1) == 1 && size(intOPs_i,2) == 4
            # If a weight vector exists and its (single) value is zero -> fallback to raw
            if !isempty(weights_i) && sum(weights_i) == 0.0
                wf.States_WF[start_idx, 1] = u[iT]
            else
                r = intOPs_i
                idx1 = Int(r[1]); w1 = r[2]
                idx2 = Int(r[3]); w2 = r[4]
                # Validate indices before accessing States_WF
                valid = 0 < idx1 <= nStatesWF && 0 < idx2 <= nStatesWF
                if valid
                    wf.States_WF[start_idx, 1] = wf.States_WF[idx1, 1] * w1 + wf.States_WF[idx2, 1] * w2
                else
                    wf.States_WF[start_idx, 1] = u[iT]  # fallback
                end
            end
        elseif size(intOPs_i,1) > 0 && size(intOPs_i,2) == 4
            if isempty(weights_i) || sum(weights_i) == 0.0
                wf.States_WF[start_idx, 1] = u[iT]
            else
                sum_wU = 0.0
                sum_w = 0.0
                nrows = size(intOPs_i,1)
                @inbounds for iiT in 1:nrows
                    w = (iiT <= length(weights_i)) ? weights_i[iiT] : 0.0
                    w == 0.0 && continue
                    row = intOPs_i[iiT, :]
                    idx1 = Int(row[1]); w1 = row[2]
                    idx2 = Int(row[3]); w2 = row[4]
                    # Validate indices before accessing States_WF
                    if !(0 < idx1 <= nStatesWF && 0 < idx2 <= nStatesWF)
                        continue
                    end
                    local_u = wf.States_WF[idx1, 1] * w1 + wf.States_WF[idx2, 1] * w2
                    sum_wU += w * local_u
                    sum_w += w
                end
                wf.States_WF[start_idx, 1] = sum_w > 0 ? (sum_wU / sum_w) : u[iT]
            end
        else
            # Fallback for malformed/missing interpolation data
            wf.States_WF[start_idx, 1] = u[iT]
        end
    end
    return wf, wind
end

"""
        getDataVel(set::Settings, wind::Wind, wf::WindFarm, t, tmp_m, floris::Floris)

Return the wind speed vector `u` for all turbines at simulation time `t` according to the
configured input / model mode. Also returns (potentially updated) `wind`.

Arguments
- `set::Settings`: simulation settings (uses `set.vel_mode`)
- `wind::Wind`: wind field state (`wind.input_vel` selects special branches)
- `wf::WindFarm`: wind farm (uses `wf.nT`, `wf.States_WF`)
- `t`: current simulation time
- `tmp_m`: temporary matrix (only used for I_and_I wake reduction)
- `floris::Floris`: FLORIS parameters (yaw exponent etc., only I_and_I branch)

Supported (unit tested in `test_getDataVel_branches.jl`)
- Standard interpolation / constant variants via `set.vel_mode`:
    `Velocity_Constant`, `Velocity_Interpolation`, `Velocity_Constant_wErrorCov`,
    `Velocity_Interpolation_wErrorCov`, `Velocity_InterpTurbine`,
    `Velocity_InterpTurbine_wErrorCov`, `Velocity_ZOH_wErrorCov`.
- EnKF turbine interpolation branch: `wind.input_vel == "EnKF_InterpTurbine"` calling
    `getWindSpeedT_EnKF(Velocity_EnKF_InterpTurbine(), ...)` with clamping of out-of-range times.

Not yet fully integrated (guarded / broken tests)
- `"I_and_I"`: an internal estimator state struct (`WSEStruct` in `windfield_velocity.jl`) already exists
    and the low-level update routine `WindSpeedEstimatorIandI_FLORIDyn` runs, but a public, documented
    construction path (export, convenience constructor, validation of required fields, tests) is missing.
    The branch is therefore kept experimental and the test remains `@test_broken` until we provide a
    stable API (e.g. `build_IandI_estimator(wf, data; kwargs...)`).
- `"RW_with_Mean"`: random-walk-with-mean model commented out; current call raises `MethodError`.

Planned cleanups / TODO
1. Provide concrete exported estimator type & finalize `I_and_I` logic.
2. Re-introduce Random Walk with Mean model (`Velocity_RW_with_Mean`).

Behavior summary
- Default: `u = getWindSpeedT(set.vel_mode, wind.vel, 1:nT, t)`.
- EnKF: per-turbine linear interpolation table with time clamping.
- I_and_I: (future) estimator integration plus optional wake reduction using `tmp_m[:,1]`.
- RW_with_Mean: (future) stochastic update around mean with mean-pull term.

All returned velocities are in m/s. Only I_and_I may mutate `wind.vel` estimator state.

Example
```julia
u, wind = getDataVel(set, wind, wf, 100.0, tmp_m, floris)
```
"""
function getDataVel(set::Settings, wind::Wind, wf::WindFarm, t, tmp_m, floris::Floris)
    idx = 1:wf.nT  # avoid temporary allocation from collect
    u = nothing
    if wind.input_vel == "I_and_I"
        u, wind.vel = getWindSpeedT(set.vel_mode, wind.vel, idx, t,
                                    wf.States_WF[wf.StartI, 2], floris.p_p)
        if (t - wind.vel.StartTime) > wind.vel.WSE.Offset
            u = u ./ tmp_m[:, 1]
        end
    elseif wind.input_vel == "RW_with_Mean"
        u = getWindSpeedT(wf.States_WF[wf.StartI, 1], wind.vel)
    elseif wind.input_vel == "EnKF_InterpTurbine"
        u = getWindSpeedT_EnKF(Velocity_EnKF_InterpTurbine(), wind.vel, idx, t)
    else
        u = getWindSpeedT(set.vel_mode, wind.vel, idx, t)
    end
    return u, wind
end
