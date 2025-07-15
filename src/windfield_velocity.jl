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


