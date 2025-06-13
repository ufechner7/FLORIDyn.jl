"""
    getWindDirT(WindDir, iT, _)

Return wind direction in SOWFA-degrees for the requested turbine(s).

# Arguments
- `WindDir`: The wind direction (scalar).
- `iT`: Index or indices of the turbines (can be an integer or array).
- `_`: Placeholder for unused argument.

# Returns
- `phi`: Array of wind direction values, same size as `iT`.
"""
function getWindDirT(WindDir, iT, _)
    return fill(WindDir, size(iT))
end

"""
    getWindDirT(WindDir, iT)

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `WindDir`: Struct with fields
    - `Data::Float64`: wind direction value
    - `CholSig::AbstractMatrix`: Cholesky factor of covariance matrix (nT x nT)
- `iT`: Vector of turbine indices (can be any indexable collection)

# Returns
- `phi`: Vector of wind directions for the selected turbines, including random perturbation
"""
function getWindDirT(WindDir, iT)
    n = length(iT)
    phi = fill(WindDir.Data, n)
    # randn(n) gives a vector of n normal random numbers
    # WindDir.CholSig should
    phi .= phi .+ WindDir.CholSig * randn(n)
    return phi
end

"""
    getWindDirT_EnKF(WindDir, iT, t)

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `WindDir::AbstractMatrix`: Matrix where each row is [time, phi_T0, phi_T1, ... phi_Tn]
- `iT`: Index or indices of the turbines (can be integer or vector)
- `t`: Time of request (scalar)

# Returns
- `phi`: Wind direction(s) at time `t` for turbine(s) `iT`
"""
function getWindDirT_EnKF(WindDir::AbstractMatrix, iT, t)
    times = WindDir[:, 1]
    n_turbines = size(WindDir, 2) - 1

    # Clamp t to within the bounds of available times
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Prepare interpolation for each turbine
    phi_out = similar(WindDir[1, 2:end])
    for j in 1:n_turbines
        itp = linear_interpolation(times, WindDir[:, j+1], extrapolation_bc=Flat())
        phi_out[j] = itp(t)
    end

    return phi_out[iT]
end

"""
    getWindDirT(WindDir, iT, t)

Returns the wind direction at the respective turbine(s).
Uniform interpolation version - all turbines experience the same changes.

Arguments:
- WindDir: Nx2 Array, columns are time and phi (wind direction)
- iT: single value or array with turbine index/indices
- t: time of request

Returns:
- phi: Array of wind directions for each turbine in iT
"""
function getWindDirT(WindDir::AbstractMatrix, iT, t)
    times = WindDir[:, 1]
    phis = WindDir[:, 2]

    # Clamp t to bounds of times
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Linear interpolation (like interp1 in MATLAB)
    itp = LinearInterpolation(times, phis, extrapolation_bc=Flat())
    phi_val = itp(t)

    # Return phi for each turbine in iT (broadcasted)
    return fill(phi_val, length(iT))
end
