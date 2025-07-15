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
