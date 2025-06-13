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
