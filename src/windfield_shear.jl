using Interpolations

"""
    getWindShearT(WindShear, z)

Compute the wind shear at a given height `z` using the specified `WindShear` model.

# Arguments
- `WindShear`: An object or parameter set describing the wind shear profile.
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
In case z IS OUT OF BOUNDS the function will use the closest available setpoint.
"""
function getWindShearT(WindShear, z)
    # Extract columns
    heights = WindShear[:, 1]
    speeds = WindShear[:, 2]

    # Handle out-of-bounds: clamp z to [minZ, maxZ]
    minZ = minimum(heights)
    maxZ = maximum(heights)
    z_clamped = clamp.(z, minZ, maxZ)

    # Linear interpolation
    itp = linear_interpolation(heights, speeds, extrapolation_bc=Flat())
    shear = itp(z_clamped)
end
