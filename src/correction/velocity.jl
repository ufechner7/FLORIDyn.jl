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

Retrieve wind velocity data for all turbines at the current simulation time using various input modes.

# Arguments
- `set::Settings`: Simulation settings containing the velocity mode configuration
- `wind::Wind`: Wind field data structure containing velocity information and input type
- `wf::WindFarm`: Wind farm object containing turbine states and configuration
- `t`: Current simulation time for temporal interpolation
- `tmp_m`: Temporary matrix containing reduction factors and other intermediate calculations
- `floris::Floris`: FLORIS model parameters including wake model coefficients

# Returns
- `u`: Wind velocity values for all turbines at the specified time
- `wind::Wind`: Updated wind field data structure (may be modified for certain input modes)

# Description
This function reads wind velocity data and returns the current wind speed (U) for all turbines 
in the wind farm. The function handles multiple input modes through conditional logic:

1. **Interpolation and Integration mode** (`"I_and_I"`): Uses `getWindSpeedT` with turbine positions, 
   wind direction, and FLORIS parameters. Applies wake effects reduction if the simulation time 
   exceeds the wind speed estimator offset.

2. **Random Walk with Mean mode** (`"RW_with_Mean"`): Uses the current wind farm state from 
   `wf.States_WF[wf.StartI, 1]` with the wind velocity model.

3. **Ensemble Kalman Filter mode** (`"EnKF_InterpTurbine"`): Uses EnKF-based interpolation 
   at turbine locations via `getWindSpeedT_EnKF`.

4. **Standard temporal interpolation mode** (default): Uses direct temporal interpolation 
   with `getWindSpeedT` for all turbines.

## Current Support Status

Supported and tested:
- Default branch with velocity model types passed via `set.vel_mode`:
    - `Velocity_Constant`
    - `Velocity_Interpolation`
    - `Velocity_Constant_wErrorCov`
    - `Velocity_Interpolation_wErrorCov`
    - `Velocity_InterpTurbine`
    - `Velocity_InterpTurbine_wErrorCov`
    - `Velocity_ZOH_wErrorCov`

Partially implemented / NOT currently operational:
- `wind.input_vel == "I_and_I"`:
    - Marked "NOT YET IMPLEMENTED" in source. Requires a dedicated estimator struct
        providing fields (`WSE`, timing, yaw, torque, pitch, etc.). Existing pathway
        executes but lacks a stable, exported type and comprehensive tests.
    - Tests are marked `@test_broken`.
- `wind.input_vel == "RW_with_Mean"`:
    - Underlying `getWindSpeedT(::Velocity_RW_with_Mean, ...)` implementation is
        commented out. Current call pattern will raise a `MethodError`.
- `wind.input_vel == "EnKF_InterpTurbine"`:
    - Dispatch for `getWindSpeedT_EnKF` expects a leading model type
        (`Velocity_EnKF_InterpTurbine`) but `getDataVel` currently calls it without
        passing that type, so the branch raises a `MethodError`.

Planned / future action items:
1. Introduce concrete estimator type(s) for I_and_I and integrate with exported API.
2. Implement / restore Random Walk with Mean velocity model (`Velocity_RW_with_Mean`).
3. Align EnKF branch call signature (e.g. `getWindSpeedT_EnKF(set.vel_mode, wind.vel, ...)`).
4. Replace ad-hoc collection allocations (`collect(1:wf.nT)`) with reused buffers.

Until these are addressed, only the default branch (and the velocity model types listed
above) should be relied upon in production runs.

# Examples
```julia
# Get wind velocity for all turbines at current simulation time
u, updated_wind = getDataVel(settings, wind_data, wind_farm, 100.0, temp_matrix, floris_params)

# The returned u contains velocity values for all turbines
velocity_turbine_1 = u[1]
```

# Notes
- The function automatically handles different wind input modes through conditional logic
- For I_and_I mode, wake effects are applied based on timing and reduction factors
- The wind structure may be modified and returned for certain input modes
- Velocity values are in m/s
- Wake reduction is only applied in I_and_I mode when sufficient time has elapsed

# See also
- [`Settings`](@ref): Simulation settings structure
- [`Wind`](@ref): Wind field data structure
- [`WindFarm`](@ref): Wind farm configuration structure
- [`Floris`](@ref): FLORIS model parameters structure
"""
function getDataVel(set::Settings, wind::Wind, wf::WindFarm, t, tmp_m, floris::Floris)
    # Initialize u
    u = nothing

    # Determine which input mode to use
    if wind.input_vel == "I_and_I"
        u, wind.vel = getWindSpeedT(set.vel_mode, wind.vel, collect(1:wf.nT), t,
                                   wf.States_WF[wf.StartI, 2], floris.p_p)

        if (t - wind.vel.StartTime) > wind.vel.WSE.Offset
            # Ufree = Ueff / reduction
            u = u ./ tmp_m[:, 1]
        end
    elseif wind.input_vel == "RW_with_Mean"
        u = getWindSpeedT(wf.States_WF[wf.StartI, 1], wind.vel)
    elseif wind.input_vel == "EnKF_InterpTurbine"
        u = getWindSpeedT_EnKF(wind.vel, collect(1:wf.nT), t)
    else
        u = getWindSpeedT(set.vel_mode, wind.vel, collect(1:wf.nT), t)
    end

    return u, wind
end

