# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Functions to calculate the wind direction depending on time and location.

# helper function
issquare(A) = size(A, 1) == size(A, 2)

"""
    getWindDirT(::Direction_Constant, WindDir, iT, _)

Return wind direction in SOWFA-degrees for the requested turbine(s).

# Arguments
- `WindDir`: The wind direction (scalar).
- `iT`: Index or indices of the turbines (can be an integer or vector).
- `_`: Placeholder for unused argument.

# Returns
- `phi`: Array of wind direction values, same size as `iT`.
"""
function getWindDirT(::Direction_Constant, WindDir, iT, _)
    if isa(iT, AbstractArray)
        return fill(WindDir, size(iT))
    else
        return WindDir
    end
end

"""
    getWindDirT(::Direction_Constant_wErrorCov, WindDir::WindDirType, iT, t)

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `WindDir::WindDirType`: [WindDirType](@ref)
- `iT`: Vector of turbine indices (can be any indexable collection)
- `t`: Time step

# Returns
- `phi`: Vector of wind directions for the selected turbines, including random perturbation
"""
function getWindDirT(::Direction_Constant_wErrorCov, WindDir::WindDirType, iT, t)
    if isa(iT, AbstractArray)
        n = length(iT)
        indices = iT
    else
        n = 1
        indices = [iT]
    end
    phi = fill(WindDir.Data, n)
    # randn(RNG,n) gives a vector of n normal random numbers
    @assert issquare(WindDir.CholSig)
    # Extract the relevant submatrix of the Cholesky factor for the selected turbines
    chol_sub = WindDir.CholSig[indices, indices]
    phi .+= chol_sub * randn(RNG, n)
    return isa(iT, AbstractArray) ? phi : phi[1]
end

# function phi = getWindDirT(WindDir,iT,~)
# %GETWINDDIRT Return wind direction in SOWFA-deg for the requested
# %turbine(s)
# %
# % ======= Input ======
# % WindDir.Data   = float, wind direction
# % WindDir.ColSig = nT x nT, col(Covariance Matrix)
# % iT        = Index/Indeces of the turbines

# phi = ones(size(iT))*WindDir.Data;
# phi = phi + (randn(1,length(phi))*WindDir.CholSig)';
# end


# # Backward compatibility method with 3 arguments
# function getWindDirT(::Direction_Constant_wErrorCov, WindDir::WindDirType, iT)
#     n = length(iT)
#     phi = fill(WindDir.Data, n)
#     # randn(RNG,n) gives a vector of n normal random numbers
#     @assert issquare(WindDir.CholSig)
#     phi .+= WindDir.CholSig * randn(RNG, n)
#     return phi
# end

"""
    getWindDirT_EnKF(::Direction_EnKF_InterpTurbine, WindDir::Matrix, iT, t)

# Direction_EnKF_InterpTurbine 

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `WindDir::Matrix`: Matrix where each row is [time, `phi_T0`, `phi_T1`, ... `phi_Tn`]
- `iT`: Index or indices of the turbines (can be integer or vector)
- `t`: Time of request (scalar)

# Returns
- `phi`: Wind direction(s) at time `t` for turbine(s) `iT`
"""
function getWindDirT_EnKF(::Direction_EnKF_InterpTurbine, WindDir::Matrix, iT, t)
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
    getWindDirT(::Direction_Interpolation, WindDir::Matrix, iT, t)

# Direction_Interpolation 

Returns the wind direction at the respective turbine(s).
Uniform interpolation version - all turbines experience the same changes.

Arguments:
- WindDir::Matrix: columns are time and phi (wind direction)
- iT: single value or vector with turbine index/indices
- t: time of request

Returns:
- phi: Vector of wind directions for each turbine in iT
"""
function getWindDirT(::Direction_Interpolation, WindDir::Matrix, iT, t)
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
    getWindDirT(::Direction_Interpolation_wErrorCov, WindDir::WindDirMatrix, iT, t)

Returns the wind direction at the respective turbine(s).
Uniform interpolation version - all turbines experience the same changes.

Arguments:
- WindDir::WindDirMatrix: [WindDirMatrix](@ref)
- iT: single value or vector with turbine index/indices
- t: time of request

Returns:
- phi: Vector of wind directions for each turbine in iT
"""
function getWindDirT(::Direction_Interpolation_wErrorCov, WindDir::WindDirMatrix, iT, t)
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

    if isa(iT, AbstractArray)
        # Array case: use full covariance matrix
        phi = fill(phi_val, length(iT))
        phi += (randn(RNG, 1, length(phi)) * WindDir.CholSig)'
        return phi
    else
        # Scalar case: use only diagonal element
        phi = phi_val
        sigma = WindDir.CholSig[iT, iT]  # Get diagonal element for this turbine
        phi += randn(RNG) * sigma
        return phi
    end
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
    getWindDirT(::Direction_InterpTurbine, WindDir::Matrix, iT, t)

Return wind direction in SOWFA-degrees for the requested turbine(s).

# Arguments
- `WindDir::Matrix`: Each row is `[time, phi_T0, phi_T1, ...]`.
- `iT: Index or indices of turbines.
- `t`: Time of request.

# Returns
- `phi::Vector{Float64}`: Wind direction(s) for the selected turbine(s) at time `t`.
"""
function getWindDirT(::Direction_InterpTurbine, WindDir, iT, t)
    # Handle case where WindDir is Nothing
    if WindDir === nothing
        error("WindDir data is missing for Direction_InterpTurbine mode. Please provide wind direction data or use a different mode.")
    end
    
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
    getWindDirT(::Direction_InterpTurbine_wErrorCov, WindDir::WindDirMatrix, iT, t)

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `WindDir::WindDirMatrix`: See: [`WindDirMatrix`](@ref)
- `iT`: Index or indices of the turbines (can be integer or vector)
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
    noise = randn(RNG,length(iT_vec))
    phi = phi .+ WindDir.CholSig[iT_vec, iT_vec] * noise

    return phi
end

"""
    getWindDirT(::Direction_RW_with_Mean, WindDirNow, WindDir::WindDirTriple)

Returns the wind direction at the respective turbine(s).

# Arguments
- `WindDirNow`: Current value (vector)
- `WindDir::WindDirTriple`: [`WindDirTriple`](@ref)

# Returns
- `phi`: Updated wind direction(s) (vector)
"""
function getWindDirT(::Direction_RW_with_Mean, WindDirNow, WindDir::WindDirTriple)
    # Random walk model with mean implementation
    # Generate random normal vector
    weightedRandN = randn(RNG,1, length(WindDirNow))
    # Compute new wind direction
    phi = WindDirNow .+ (weightedRandN * WindDir.CholSig)' .+
          WindDir.MeanPull .* (WindDir.Init .- WindDirNow)
    return phi
end

"""
    getWindDirT(::Direction_RW_with_Mean, WindDir::WindDirTriple, iT, t)

Random walk with mean reversion model for wind direction.

# Arguments
- `::Direction_RW_with_Mean`: Direction mode indicator
- `WindDir::WindDirTriple`: Wind direction data containing Init, CholSig, and MeanPull
- `iT`: Turbine index or indices
- `t`: Time value (unused in this implementation)

# Returns
- `phi`: Wind direction(s) for the requested turbine(s)
"""
function getWindDirT(::Direction_RW_with_Mean, WindDir::WindDirTriple, iT, t)
    if isa(iT, AbstractArray)
        n = length(iT)
        indices = iT
    else
        n = 1
        indices = [iT]
    end
    
    # Get initial values for selected turbines
    WindDirNow = WindDir.Init[indices]
    
    # Generate random normal vector
    weightedRandN = randn(RNG, 1, n)
    
    # Extract relevant submatrix for selected turbines
    chol_sub = WindDir.CholSig[indices, indices]
    
    # Compute new wind direction with mean reversion
    phi = WindDirNow .+ (weightedRandN * chol_sub)' .+
          WindDir.MeanPull .* (WindDir.Init[indices] .- WindDirNow)
    
    return isa(iT, AbstractArray) ? phi : phi[1]
end
