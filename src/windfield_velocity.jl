# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getWindSpeedT(::Velocity_Constant, WindVel::Number, iT)

Returns the wind speed at a given time index `iT` for a constant velocity wind field.

# Arguments
- `::Velocity_Constant`: Type indicator for constant wind velocity.
- `WindVel::Number`: The scalar wind velocity.
- `iT`: Single value or array with turbine index/indices.

# Returns
- The wind speed for the respective turbine(s), scalar or vector.
"""
function getWindSpeedT(::Velocity_Constant, WindVel, iT)
    WindVel .* ones(eltype(WindVel), size(iT))
end

"""
    getWindSpeedT(::Velocity_Constant_wErrorCov, WindVel::WindVelType, iT, _)

Compute the wind speed at a given time step `iT` using a constant velocity model with error covariance.

# Arguments
- `::Velocity_Constant_wErrorCov`: Type indicating the constant velocity model with error covariance.
- `WindVel::WindVelType`: [WindVelType](@ref)
- `iT`: Scalar or array of turbine indices
- `_`: Placeholder for additional arguments (unused).

# Returns
- The wind speed at the respective turbine(s)
"""
function getWindSpeedT(::Velocity_Constant_wErrorCov, WindVel::WindVelType, iT)
    # Create a vector of the same size as iT filled with WindVel.Data
    Vel = fill(WindVel.Data, length(iT))

    # Add Gaussian noise scaled by WindVel.CholSig
    Vel .+= (randn(1, length(Vel)) * WindVel.CholSig)'

    return Vel
end

using Interpolations

"""
    getWindSpeedT_EnKF(::Velocity_EnKF_InterpTurbine, WindVel::Matrix, iT, t)

Returns the wind speed at the turbine index or indices `iT` at time `t`, using linear interpolation of the measurement table `WindVel`. 

- WindVel: A matrix of size (ntimes, nturbines+1). First column is time, rest are wind speeds for each turbine.
- iT: Index or array of indices of turbines for output.
- t: Time at which to query wind speed.
"""
function getWindSpeedT_EnKF(::Velocity_EnKF_InterpTurbine, WindVel::Matrix, iT, t)
    times  = WindVel[:, 1]
    speeds = WindVel[:, 2:end]  # size: (N_time, N_turbines)

    # Clamp time to in-bounds with warning
    mint, maxt = times[1], times[end]
    if t < mint
        @warn "The time $t is out of bounds, will use $mint instead."
        t = mint
    elseif t > maxt
        @warn "The time $t is out of bounds, will use $maxt instead."
        t = maxt
    end

    # Interpolate each turbine column independently
    wind_at_t = [LinearInterpolation(times, speeds[:, j], extrapolation_bc=Flat())(t) for j in 1:size(speeds, 2)]

    # Select desired turbine(s)
    return wind_at_t[iT]
end



