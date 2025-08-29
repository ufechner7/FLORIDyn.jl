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

Influence-based free-stream wind speed correction translated from the legacy MATLAB
`correctVel` routine. Each turbine's ambient (free) wind speed state `States_WF[:,1]`
can depend on upstream observation point (OP) interpolations defined in `wf.intOPs`
and weighted by wake influence factors in `wf.Weight` with dependency lists in `wf.dep`.

# Behavior per turbine iT
1. No dependencies (`wf.dep[iT]` empty): assign raw velocity `u[iT]`.
2. Single interpolation row (size 1x4): treat row as `[idx1 w1 idx2 w2]` and set
   velocity to `w1*U[idx1] + w2*U[idx2]` using already updated `wf.States_WF[idx,1]`.
3. Multiple rows (each 4 columns) with non-zero total weight: compute per-row
   interpolated velocities then take weighted average using `wf.Weight[iT]`.
4. Missing / zero weights -> fallback to raw `u[iT]`.
5. Malformed or missing interpolation data -> fallback to raw `u[iT]`.

# Arguments
- `::Velocity_Influence`: Strategy marker dispatch type
- `set::Settings`: Simulation settings (used to obtain base velocities via `getDataVel`)
- `wf::WindFarm`: Wind farm state (modified in-place; column 1 is free wind speed)
- `wind::Wind`: Wind data source (may be updated by `getDataVel` in estimator modes)
- `t`: Simulation time
- `floris`: FLORIS parameter container (passed through to `getDataVel` if needed)
- `tmpM`: Temporary matrix used in special velocity retrieval modes (e.g. I_and_I)

# Returns
- Updated `(wf, wind)` tuple (mutates `wf` in-place; returns `wind` for API symmetry)

# Notes
- Assumes rows of `wf.intOPs[iT]` are `[idx1 w1 idx2 w2]`.
- Uses already updated values within the same loop for chained dependencies.
- Robust to missing `dep`, `Weight` or `intOPs` entries (graceful fallback).
"""
function correctVel(::Velocity_Influence, set::Settings, wf::WindFarm, wind::Wind, t, floris, tmpM)
    # Base free wind speeds (may update wind state depending on mode)
    u, wind = getDataVel(set, wind, wf, t, tmpM, floris)

    nT = wf.nT
    has_dep = !isempty(wf.dep)
    has_intOPs = !isempty(wf.intOPs)
    has_weights = !isempty(wf.Weight)

    for iT in 1:nT
        dep_i = (has_dep && length(wf.dep) >= iT) ? wf.dep[iT] : Int[]
        start_idx = wf.StartI[iT, 1]

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
                wf.States_WF[start_idx, 1] = wf.States_WF[idx1, 1] * w1 + wf.States_WF[idx2, 1] * w2
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

