# the following paths are used:
# - Direction_EnKF_InterpTurbine



"""
    getWindDirT(::Direction_Constant, WindDir, iT, _)

Return wind direction in SOWFA-degrees for the requested turbine(s).

# Arguments
- `WindDir`: The wind direction (scalar).
- `iT`: Index or indices of the turbines (can be an integer or array).
- `_`: Placeholder for unused argument.

# Returns
- `phi`: Array of wind direction values, same size as `iT`.
"""
function getWindDirT(::Direction_Constant, WindDir, iT, _)
    return fill(WindDir, size(iT))
end

"""
    getWindDirT(::Direction_Constant_wErrCov, WindDir, iT)

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `WindDir`: Struct with fields
    - `Data::Float64`: wind direction value
    - `CholSig::AbstractMatrix`: Cholesky factor of covariance matrix (nT x nT)
- `iT`: Vector of turbine indices (can be any indexable collection)

# Returns
- `phi`: Vector of wind directions for the selected turbines, including random perturbation
"""
function getWindDirT(::Direction_Constant_wErrCov, WindDir, iT)
    n = length(iT)
    phi = fill(WindDir.Data, n)
    # randn(n) gives a vector of n normal random numbers
    # WindDir.CholSig should
    phi .= phi .+ WindDir.CholSig * randn(n)
    return phi
end

"""
    getWindDirT_EnKF(::Direction_EnKF_InterpTurbine, WindDir, iT, t)

# Direction_EnKF_InterpTurbine 

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `WindDir::AbstractMatrix`: Matrix where each row is [time, phi_T0, phi_T1, ... phi_Tn]
- `iT`: Index or indices of the turbines (can be integer or vector)
- `t`: Time of request (scalar)

# Returns
- `phi`: Wind direction(s) at time `t` for turbine(s) `iT`
"""
function getWindDirT_EnKF(::Direction_EnKF_InterpTurbine, WindDir::AbstractMatrix, iT, t)
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
    getWindDirT(::Direction_Interpolation, WindDir, iT, t)

# Direction_Interpolation 

Returns the wind direction at the respective turbine(s).
Uniform interpolation version - all turbines experience the same changes.

Arguments:
- WindDir: Nx2 Array, columns are time and phi (wind direction)
- iT: single value or array with turbine index/indices
- t: time of request

Returns:
- phi: Array of wind directions for each turbine in iT
"""
function getWindDirT(::Direction_Interpolation, WindDir::AbstractMatrix, iT, t)
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
    itp = linear_interpolation(times, phis, extrapolation_bc=Flat())
    phi_val = itp(t)

    # Return phi for each turbine in iT (broadcasted)
    return fill(phi_val, length(iT))
end

"""
    getWindDirT(::Direction_Interpolation_wErrorCov, WindDir, iT, t)

Returns the wind direction at the respective turbine(s).
Uniform interpolation version - all turbines experience the same changes.

Arguments:
- WindDir.Data   = (t,phi) pairs between which is linearly interpolated
- WindDir.ColSig = nT x nT, col(Covariance Matrix)
- iT: single value or array with turbine index/indices
- t: time of request

Returns:
- phi: Array of wind directions for each turbine in iT
"""
function getWindDirT(::Direction_Interpolation_wErrorCov, WindDir, iT, t)
    # Ensure t is within bounds
    if t < WindDir.Data[1, 1]
        @warn "The time $t is out of bounds, will use $(WindDir.Data[1,1]) instead."
        t = WindDir.Data[1, 1]
    elseif t > WindDir.Data[end, 1]
        @warn "The time $t is out of bounds, will use $(WindDir.Data[end,1]) instead."
        t = WindDir.Data[end, 1]
    end

    # Linear interpolation (equivalent to interp1 in MATLAB)
    phi_val = interp(WindDir.Data[:, 1], WindDir.Data[:, 2], t)

    # Replicate for all turbines in iT
    phi = fill(phi_val, length(iT))

    # Add correlated noise (assuming CholSig is a Cholesky factor)
    phi += (randn(1, length(phi)) * WindDir.CholSig)'

    return phi
end

# Helper function for 1D linear interpolation
function interp(x, y, t)
    idx = searchsortedlast(x, t)
    if idx == length(x)
        return y[end]
    elseif idx == 0
        return y[1]
    else
        x0, x1 = x[idx], x[idx+1]
        y0, y1 = y[idx], y[idx+1]
        return y0 + (y1 - y0) * (t - x0) / (x1 - x0)
    end
end

"""
    getWindDirT(::Direction_InterpTurbine, WindDir, iT, t)

Return wind direction in SOWFA-degrees for the requested turbine(s).

# Arguments
- `WindDir::Matrix{Float64}`: Each row is `[time, phi_T0, phi_T1, ...]`.
- `iT::Union{Int, AbstractVector{Int}}`: Index or indices of turbines (1-based).
- `t::Float64`: Time of request.

# Returns
- `phi::Vector{Float64}`: Wind direction(s) for the selected turbine(s) at time `t`.
"""
function getWindDirT(::Direction_InterpTurbine, WindDir::Matrix{Float64}, iT, t::Float64)
    # Check time bounds
    tmin = WindDir[1, 1]
    tmax = WindDir[end, 1]
    if t < tmin
        @warn "The time $t is out of bounds, will use $tmin instead."
        t = tmin
    elseif t > tmax
        @warn "The time $t is out of bounds, will use $tmax instead."
        t = tmax
    end

    # Interpolate for all turbines at time t
    times = WindDir[:, 1]
    phis = WindDir[:, 2:end]  # Each column is a turbine
    phi_out = [interp(times, phis[:, j], t) for j in 1:size(phis, 2)]

    # Select requested turbines
    return phi_out[iT]
end

"""
    getWindDirT(::Direction_InterpTurbine_wErrorCov, WindDir, iT, t)

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `WindDir`: A struct with fields:
    - `Data` (Matrix): Each row is [time, phi_T0, phi_T1, ... phi_Tn]
    - `CholSig` (Matrix): Cholesky factor of the covariance matrix
- `iT`: Index or indices of the turbines (can be integer or array)
- `t`: Time of request (Float64)

# Returns
- `phi`: Wind direction(s) for requested turbine(s), perturbed with noise.
"""
function getWindDirT(::Direction_InterpTurbine_wErrorCov, WindDir, iT, t)
    times = WindDir.Data[:, 1]
    nTurbines = size(WindDir.Data, 2) - 1

    # Clamp t to available time range, with warning
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Interpolate for each turbine at time t
    phi_out = [Interpolations.linear_interpolation(times, WindDir.Data[:, j+1])(t) for j in 1:nTurbines]

    # Select requested turbines
    phi = phi_out[iT]

    # Add correlated noise
    iT_vec = isa(iT, Integer) ? [iT] : iT
    phi = phi_out[iT_vec]
    noise = randn(length(iT_vec))
    phi = phi .+ WindDir.CholSig[iT_vec, iT_vec] * noise

    return phi
end


