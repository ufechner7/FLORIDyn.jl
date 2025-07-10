# function shear = getWindShearT(WindShear,z)
# %GETWINDSHEART returns the relative reduction (or speed-up) by wind shear
# % Expects a .csv file called "WindShearProfile.csv" with a normalized wind
# % speed profile for different heights:
# %   z, (u_z/u0)
# %   z, (u_z/u0)
# %   z, (u_z/u0)
# % There is a linear interpolation between every pair
# % IN CASE z IS OUT OF BOUNDS the function will use the closest available
# % setpoint
# % ======================================================================= %
# % WindShear = normalized wind speed at different heights
# % z         = height(s)
# % ======================================================================= %
# % Out of bounds handling
# maxZ = max(WindShear(:,1));
# minZ = min(WindShear(:,1));
# z(z>maxZ) = maxZ;
# z(z<minZ) = minZ;
# % Interpolate
# shear = interp1(WindShear(:,1),WindShear(:,2),z);
# end

using Interpolations

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

    return shear
end
