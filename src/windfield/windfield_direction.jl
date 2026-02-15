# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Functions to calculate the wind direction depending on time and location.

# helper function
issquare(A) = size(A, 1) == size(A, 2)

"""
    getWindDirT(::Direction_Constant, wind::Wind, iT, _)

Return wind direction in SOWFA-degrees for the requested turbine(s).

# Arguments
- `wind`: The wind direction (scalar).
- `iT`: Index or indices of the turbines (can be an integer or vector).
- `_`: Placeholder for unused argument.

# Returns
- `phi`: Array of wind direction values, same size as `iT`.
"""
function getWindDirT(::Direction_Constant, wind::Wind, iT, _)
    wind_dir = wind.dir_fixed
    if isa(iT, AbstractArray)
        return fill(wind_dir, size(iT))
    else
        return wind_dir
    end
end

"""
    getWindDirT(::Direction_Constant_wErrorCov, wind::Wind, iT, t)

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `wind::Wind`: [Wind](@ref)
- `iT`: Vector of turbine indices (can be any indexable collection)
- `t`: Time step

# Returns
- `phi`: Vector of wind directions for the selected turbines, including random perturbation
"""
function getWindDirT(::Direction_Constant_wErrorCov, wind::Wind, iT, t)
    wind_dir = wind.dir
    if isa(iT, AbstractArray)
        n = length(iT)
        indices = iT
    else
        n = 1
        indices = [iT]
    end
    phi = fill(wind_dir.Data, n)
    # randn(RNG,n) gives a vector of n normal random numbers
    @assert issquare(wind_dir.CholSig)
    # Extract the relevant submatrix of the Cholesky factor for the selected turbines
    chol_sub = wind_dir.CholSig[indices, indices]
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
# % iT        = Index/Indices of the turbines

# phi = ones(size(iT))*WindDir.Data;
# phi = phi + (randn(1,length(phi))*WindDir.CholSig)';
# end

"""
    getWindDirT_EnKF(::Direction_EnKF_InterpTurbine, wind, iT, t)

# Direction_EnKF_InterpTurbine 

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `wind_dir::Matrix`: Matrix where each row is [time, `phi_T0`, `phi_T1`, ... `phi_Tn`]
- `iT`: Index or indices of the turbines (can be integer or vector)
- `t`: Time of request (scalar)

# Returns
- `phi`: Wind direction(s) at time `t` for turbine(s) `iT` [°]
"""
function getWindDirT_EnKF(::Direction_EnKF_InterpTurbine, wind::Wind, iT, t)
    wind_dir = wind.dir
    times = wind_dir[:, 1]
    n_turbines = size(wind_dir, 2) - 1

    # Clamp t to within the bounds of available times
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Prepare interpolation for each turbine
    phi_out = similar(wind_dir[1, 2:end])
    for j in 1:n_turbines
        itp = linear_interpolation(times, wind_dir[:, j+1], extrapolation_bc=Flat())
        phi_out[j] = itp(t)
    end

    return phi_out[iT]
end

"""
    getWindDirT(::Direction_Interpolation, wind::Wind, iT, t)

# Direction_Interpolation 

Returns the wind direction at the respective turbine(s).
Uniform interpolation version - all turbines experience the same changes.

Arguments:
- wind_dir::Matrix: columns are time and phi (wind direction)
- iT: single value or vector with turbine index/indices
- t: time of request

Returns:
- phi: Vector of wind directions for each turbine in iT [°]
"""
function getWindDirT(::Direction_Interpolation, wind::Wind, iT, t)
    wind_dir = wind.dir
    times = wind_dir[:, 1]
    phis = wind_dir[:, 2]

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
    getWindDirT(::Direction_Interpolation_wErrorCov, wind::Wind, iT, t)

Returns the wind direction at the respective turbine(s).
Uniform interpolation version - all turbines experience the same changes.

Arguments:
- wind::Wind: [Wind](@ref)
- iT: single value or vector with turbine index/indices
- t: time of request

Returns:
- phi: Vector of wind directions for each turbine in iT [°]
"""
function getWindDirT(::Direction_Interpolation_wErrorCov, wind::Wind, iT, t)
    wind_dir = wind.dir

    # Ensure t is within bounds
    if t < wind_dir.Data[1, 1]
        @warn "The time $t is out of bounds, will use $(wind_dir.Data[1,1]) instead."
        t = wind_dir.Data[1, 1]
    elseif t > wind_dir.Data[end, 1]
        @warn "The time $t is out of bounds, will use $(wind_dir.Data[end,1]) instead."
        t = wind_dir.Data[end, 1]
    end

    # Linear interpolation (equivalent to interp1 in MATLAB)
    phi_val = interp(wind_dir.Data[:, 1], wind_dir.Data[:, 2], t)

    if isa(iT, AbstractArray)
        # Array case: use full covariance matrix
        phi = fill(phi_val, length(iT))
        phi += (randn(RNG, 1, length(phi)) * wind_dir.CholSig)'
        return phi
    else
        # Scalar case: use only diagonal element
        phi = phi_val
        sigma = wind_dir.CholSig[iT, iT]  # Get diagonal element for this turbine
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
    getWindDirT(::Direction_InterpTurbine, wind::Wind, iT, t)

Return wind direction in SOWFA-degrees for the requested turbine(s).

# Arguments
- `wind_dir::Matrix`: Each row is `[time, phi_T0, phi_T1, ...]`.
- `iT: Index or indices of turbines.
- `t`: Time of request. [s]

# Returns
- `phi::Vector{Float64}`: Wind direction(s) for the selected turbine(s) at time `t`. [°]
"""
function getWindDirT(::Direction_InterpTurbine, wind::Wind, iT, t)
    wind_dir = wind.dir

    # Handle case where wind_dir is Nothing
    if wind_dir === nothing
        error("wind_dir data is missing for Direction_InterpTurbine mode. Please provide wind direction data or use a different mode.")
    end
    
    # Check time bounds
    tmin = wind_dir[1, 1]
    tmax = wind_dir[end, 1]
    if t < tmin
        @warn "The time $t is out of bounds, will use $tmin instead."
        t = tmin
    elseif t > tmax
        @warn "The time $t is out of bounds, will use $tmax instead."
        t = tmax
    end

    # Interpolate for all turbines at time t
    times = wind_dir[:, 1]
    phis = wind_dir[:, 2:end]  # Each column is a turbine
    phi_out = [interp(times, phis[:, j], t) for j in axes(phis, 2)]

    # Select requested turbines
    return phi_out[iT]
end

"""
    getWindDirT(::Direction_InterpTurbine_wErrorCov, wind::Wind, iT, t)

Return wind direction in SOWFA-deg for the requested turbine(s).

# Arguments
- `wind::Wind`: See: [`Wind`](@ref)
- `iT`: Index or indices of the turbines (can be integer or vector)
- `t`: Time of request (Float64) [s]

# Returns
- `phi`: Wind direction(s) for requested turbine(s), perturbed with noise. [°]
"""
function getWindDirT(::Direction_InterpTurbine_wErrorCov, wind::Wind, iT, t)
    wind_dir = wind.dir
    times = wind_dir.Data[:, 1]
    nTurbines = size(wind_dir.Data, 2) - 1

    # Clamp t to available time range, with warning
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Interpolate for each turbine at time t
    phi_out = [Interpolations.linear_interpolation(times, wind_dir.Data[:, j+1])(t) for j in 1:nTurbines]

    # Select requested turbines
    phi = phi_out[iT]

    # Add correlated noise
    iT_vec = isa(iT, Integer) ? [iT] : iT
    phi = phi_out[iT_vec]
    noise = randn(RNG,length(iT_vec))
    phi = phi .+ wind_dir.CholSig[iT_vec, iT_vec] * noise

    return phi
end

"""
    getWindDirT(::Direction_RW_with_Mean, wind_dir_now, wind::Wind)

Returns the wind direction at the respective turbine(s).

# Arguments
- `wind_dir_now`: Current value (vector)
- `wind_dir::Wind`: [`Wind`](@ref)

# Returns
- `phi`: Updated wind direction(s) (vector) [°]
"""
function getWindDirT(::Direction_RW_with_Mean, wind_dir_now, wind_dir::Wind)
    wind_dir_triple = wind_dir.dir    
    # Random walk model with mean implementation
    # Generate random normal vector
    weightedRandN = randn(RNG,1, length(wind_dir_now))
    # Compute new wind direction
    phi = wind_dir_now .+ (weightedRandN * wind_dir_triple.CholSig)' .+
          wind_dir_triple.MeanPull .* (wind_dir_triple.Init .- wind_dir_now)
    return phi
end

"""
    getWindDirT(::Direction_RW_with_Mean, wind::Wind, iT, t)

Calculate wind direction using a random walk with mean reversion model.

This model combines stochastic variability with a tendency to return to a long-term 
average direction, making it suitable for modeling realistic wind direction behavior 
that exhibits both short-term fluctuations and long-term meteorological patterns.

# Model Description
The wind direction is updated according to:
```
φ(t+1) = φ(t) + random_perturbation + mean_pull × (φ_target - φ(t))
```

Where:
- `random_perturbation`: Gaussian noise with covariance structure
- `mean_pull`: Strength of reversion toward the target direction (0 = no reversion, 1 = strong reversion)
- `φ_target`: Long-term average or equilibrium wind direction

# Arguments
- `::Direction_RW_with_Mean`: Direction mode indicator
- `wind::Wind`: Wind configuration containing direction data as [`WindDirTriple`](@ref)
  - `wind.dir.Init`: Target/equilibrium wind direction(s) [°]
  - `wind.dir.CholSig`: Cholesky factor of covariance matrix for random perturbations
  - `wind.dir.MeanPull`: Mean reversion strength parameter (0-1)
- `iT`: Turbine index or indices (Integer or AbstractArray)
- `t`: Time value [s]. Note: Not used in this memoryless implementation.

# Returns
- `phi`: Wind direction(s) for the requested turbine(s) [°]

# Physical Interpretation
- **Mean reversion**: Models the tendency of wind to return to prevailing directions
- **Stochastic component**: Captures short-term meteorological variability
- **Spatial correlation**: Through the covariance matrix, accounts for spatial dependencies between turbines

# Example
```julia
# Wind tends to revert to westerly (270°) with moderate strength
wind_dir_triple = WindDirTriple(
    Init=[270.0, 270.0, 270.0],     # Target directions for 3 turbines
    MeanPull=0.1,                   # 10% reversion strength per time step
    CholSig=0.5*I(3)                # Independent 0.5° standard deviation
)
```
"""
function getWindDirT(::Direction_RW_with_Mean, wind::Wind, iT, t)
    wind_dir_triple = wind.dir
    if isa(iT, AbstractArray)
        n = length(iT)
        indices = iT
    else
        n = 1
        indices = [iT]
    end
    
    # Get initial values for selected turbines
    wind_dir_now = wind_dir_triple.Init[indices]

    # Generate random normal vector
    weightedRandN = randn(RNG, 1, n)
    
    # Extract relevant submatrix for selected turbines
    chol_sub = wind_dir_triple.CholSig[indices, indices]
    
    # Compute new wind direction with mean reversion
    phi = wind_dir_now .+ (weightedRandN * chol_sub)' .+
          wind_dir_triple.MeanPull .* (wind_dir_triple.Init[indices] .- wind_dir_now)
    
    return isa(iT, AbstractArray) ? phi : phi[1]
end
