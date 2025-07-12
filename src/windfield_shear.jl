# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Interpolations

"""
    getWindShearT(::Shear_Interpolation, WindShear::AbstractMatrix, z)

Compute the wind shear at a given height `z` using the specified `WindShear` model.

# Arguments
- `::Shear_Interpolation`: (Type only) Use interpolation to determine the wind shear.
- `WindShear`: A matrix describing the wind shear profile.
- `z`: The height (in meters) at which to evaluate the wind shear.

# Returns
- The wind shear value at height `z`.

# REMARKS
Expects a .csv file called "WindShearProfile.csv" with a normalized wind speed profile for different heights:
```
z, (u_z/u0)
z, (u_z/u0)
z, (u_z/u0)
```
There is a linear interpolation between every pair.
In case z is out of bounds the function will use the closest available setpoint.
"""
function getWindShearT(::Shear_Interpolation, wind_shear::AbstractMatrix, z)
    # Extract columns
    heights = wind_shear[:, 1]
    speeds  = wind_shear[:, 2]

    # Handle out-of-bounds: clamp z to [minZ, maxZ]
    minZ = minimum(heights)
    maxZ = maximum(heights)
    z_clamped = clamp.(z, minZ, maxZ)

    # Linear interpolation
    itp = linear_interpolation(heights, speeds, extrapolation_bc=Flat())
    shear = itp(z_clamped)
end

# function shear = getWindShearT(WindShear,z_norm)
# %GETWINDSHEART Return the shear factor u_eff = shear * u_referenceHight
# % POWER LAW implementation
# %   expects a WindShearPowerLaw.csv with the shear coefficient
# % 
# % ======================================================================= %
# % WindShear = Holds shear coefficient and reference height
# %          .z0      = reference height
# %          .alpha   = shear coefficient
# % z         = height(s)
# % ======================================================================= %
# shear = (z_norm).^WindShear.alpha;
# end

"""
    getWindShearT(::Shear_PowerLaw, WindShear, z_norm)

Return the shear factor `u_eff = shear * u_referenceHeight` using the power law.

# Arguments
- `Shear_PowerLaw`: (type only, unused) Specifies that this method applies to the power law model
- `WindShear`: A struct of type (`WindShear`)(@ref)
    - `z0`: Reference height (not used in this function)
    - `alpha`: Shear coefficient
- `z_norm`: Height(s) (can be scalar or array)

# Returns
- `shear`: The shear factor at the given height(s)
"""
function getWindShearT(::Shear_PowerLaw, wind_shear::WindShear, z_norm)
    shear = z_norm .^ wind_shear.alpha
    return shear
end
