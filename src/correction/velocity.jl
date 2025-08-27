# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function correctVel(::Velocity_None, set, wf, Wind, SimTime, paramFLORIS, tmpM)
    # Get data
    U, Wind = getDataVel(set, Wind, wf, SimTime, tmpM, paramFLORIS)

    # Correct Velocity
   wf.States_WF[wf.StartI, 1] = U

    return wf, Wind
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

Not yet implemented (guarded / broken tests)
- `"I_and_I"`: requires a stable estimator struct (state, offsets, yaw, torque, pitch). Current path
    is placeholder; test marked `@test_broken`.
- `"RW_with_Mean"`: random-walk-with-mean model commented out; current call raises `MethodError`.

Planned cleanups / TODO
1. Provide concrete exported estimator type & finalize I_and_I logic.
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

