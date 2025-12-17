# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Functions to calculate the wind velocity depending on time and location.

"""
    getWindSpeedT_EnKF(::Velocity_EnKF_InterpTurbine, wind_vel::AbstractMatrix, iT, t)

Returns the wind speed at the turbine index or indices `iT` at time `t`, using linear interpolation of 
the measurement table `wind_vel`. 

- wind_vel: A matrix of size (ntimes, nturbines+1). First column is time, rest are wind speeds for each turbine.
- iT: Index or array of indices of turbines for output.
- t: Time at which to query wind speed.
"""
function getWindSpeedT_EnKF(::Velocity_EnKF_InterpTurbine, wind_vel::AbstractMatrix, iT, t)
    times = wind_vel[:, 1]
    n_turbines = size(wind_vel, 2) - 1

    # Clamp t to within the bounds of available times
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Prepare interpolation for each turbine
    U_out = similar(wind_vel[1, 2:end])
    for j in 1:n_turbines
        itp = linear_interpolation(times, wind_vel[:, j+1], extrapolation_bc=Flat())
        U_out[j] = itp(t)
    end

    return U_out[iT]
end

"""
    getWindSpeedT(::Velocity_Constant, wind_vel::Number, iT, _)

Returns the wind speed at a given time index `iT` for a constant velocity wind field.

# Arguments
- `::Velocity_Constant`: Type indicator for constant wind velocity.
- `wind_vel::Number`: The scalar wind velocity.
- `iT`: Single value or array with turbine index/indices.
- `_`: unused

# Returns
- The wind speed for the respective turbine(s), scalar or vector.
"""
function getWindSpeedT(::Velocity_Constant, wind_vel, iT, _)
    if isa(iT, AbstractArray)
        return wind_vel .* ones(eltype(wind_vel), size(iT))
    else
        return wind_vel
    end
end

"""
    getWindSpeedT(::Velocity_Constant_wErrorCov, wind_vel::WindVelType, iT, _)

Compute the wind speed at a given time step `iT` using a constant velocity model with error covariance.

# Arguments
- `::Velocity_Constant_wErrorCov`: Type indicating the constant velocity model with error covariance.
- `wind_vel::WindVelType`: [WindVelType](@ref)
- `iT`: Scalar or array of turbine indices
- `_`: Placeholder for additional arguments (unused).

# Returns
- The wind speed at the respective turbine(s)
"""
function getWindSpeedT(::Velocity_Constant_wErrorCov, wind_vel::WindVelType, iT)
    # Create a vector of the same size as iT filled with wind_vel.Data
    Vel = fill(wind_vel.Data, length(iT))

    # Add Gaussian noise scaled by wind_vel.CholSig
    Vel .+= (randn(RNG, 1, length(Vel)) * wind_vel.CholSig)'

    return Vel
end

mutable struct TurbineProps
    FluidDensity::Float64
    RotorRadius::Float64
    GearboxRatio::Float64
    GearboxEff::Float64
    InertiaTotal::Float64
    CpFun::Function
end

mutable struct WSEStruct
    V::Vector{Float64}
    Ee::Vector{Float64}
    omega::Vector{Float64}
    beta::Float64
    gamma::Float64
    dt_SOWF::Float64
    T_prop::TurbineProps
end

function WindSpeedEstimatorIandI_FLORIDyn(WSE, Rotor_Speed, Blade_pitch, Gen_Torque, yaw, p_p)
    # === Wind Speed Estimator (I&I Method) ===
    # Based on Liu et al., IEEE Control Systems Letters, 2021
    #
    # Inputs:
    #   WSE: Wind speed estimator state (mutable object)
    #   Rotor_Speed: Rotor speed [rpm]
    #   Blade_pitch: Blade pitch [deg]
    #   Gen_Torque: Generator torque [Nm]
    #   yaw: Yaw angle difference [rad]
    #   p_p: Exponent for yaw correction
    #
    # Returns:
    #   V_out: Estimated wind speed
    #   WSE: Updated estimator state

    # === Unpack turbine and estimator parameters ===
    ρ = WSE.T_prop.FluidDensity          # Air density [kg/m^3]
    R = WSE.T_prop.RotorRadius           # Rotor radius [m]
    A = π * R^2                          # Rotor swept area [m^2]
    γ = WSE.gamma                        # Estimator gain γ
    β = WSE.beta                         # Estimator gain β
    G = WSE.T_prop.GearboxRatio          # Gearbox ratio
    I = WSE.T_prop.InertiaTotal          # Total inertia
    η = WSE.T_prop.GearboxEff            # Gearbox efficiency
    dt = WSE.dt_SOWF                     # Simulation time step [s]

    # === Calculate Tip Speed Ratio ===
    TSR = (Rotor_Speed * π * R) ./ (WSE.V * 30)  # Note: Rotor_Speed in RPM

    # === Evaluate power coefficient and correct for yaw ===
    Cp_raw = WSE.T_prop.CpFun(TSR, Blade_pitch) # Assumes interpolation func: CpFun(tsr, pitch)
    Cp = max.(Cp_raw, 0.0) .* ρ .* cos.(yaw).^p_p

    # === Handle NaN in Cp ===
    if any(isnan.(Cp))
        println("Warning: Cp out of range. TSR = $TSR, Pitch = $Blade_pitch deg.")
        println("Assuming wind speed estimate remains unchanged.")
    else
        # === Estimate aerodynamic torque ===
        aerodynamicTorque = 0.5 * ρ * A .* (WSE.V.^3) ./ (Rotor_Speed * π / 30) .* Cp

        # Clip torque to be non-negative
        aerodynamicTorque = max.(aerodynamicTorque, 0.0)

        # === Update estimator state (Extended I&I) ===
        omega_dot = (-(η * Gen_Torque .* G .- aerodynamicTorque)) ./ I
        WSE.omega .= WSE.omega .+ dt .* omega_dot
        diff_omega = -WSE.omega .+ Rotor_Speed .* π ./ 30
        WSE.Ee     .= WSE.Ee .+ dt .* diff_omega
        WSE.V      .= β .* WSE.Ee .+ γ .* diff_omega
    end

    return WSE.V, WSE
end

function getWindSpeedT(::Velocity_I_and_I, wind_vel, iT, t, wind_dir, p_p)
    # Returns the wind speed at the respective turbine(s)
    # iT        = single value or array with turbine index/indices
    # wind_vel   = Data for I&I wind speed estimator
    # t   = Current simulation time
    # Returns:
    # U         = Velocity
    # wind_vel   = Updated struct for I&I estimator
    
    if t == wind_vel.StartTime
        U = wind_vel.WSE.V[iT]
        return U, wind_vel
    end

    # Calculate current and previous index
    indCurr = round(Int, (t - wind_vel.StartTime) / wind_vel.WSE.dt_SOWF)
    indPrev = Int((wind_vel.TimePrev - wind_vel.StartTime) / wind_vel.WSE.dt_SOWF + 1)

    # Number of Turbines
    nT = wind_vel.WSE.nT

    # Run I&I estimator from last time step to current one
    for i in indPrev:indCurr
        try
            idx = (1:nT) .+ nT * (i-1)
            Rotor_Speed = wind_vel.WSE.rotorSpeed[idx, 3]
        catch
            println("Failed to get Wind Speed")
            # Optionally break or continue, depending on the desired error handling
        end
        Blade_pitch = wind_vel.WSE.bladePitch[idx, 3]
        Gen_Torque  = wind_vel.WSE.genTorque[idx, 3]

        yaw = deg2rad.(wind_dir .- wind_vel.WSE.nacelleYaw[idx, 3])

        # Run estimator (assumes WindSpeedEstimatorIandI_FLORIDyn is implemented elsewhere)
        V_out, wind_vel.WSE = WindSpeedEstimatorIandI_FLORIDyn(
            wind_vel.WSE, Rotor_Speed, Blade_pitch, Gen_Torque, yaw, p_p
        )
    end

    # Update previous time
    wind_vel.TimePrev = indCurr * wind_vel.WSE.dt_SOWF + wind_vel.StartTime

    if (t - wind_vel.StartTime) > wind_vel.WSE.Offset
        U = wind_vel.WSE.V[iT]
    else
        U = wind_vel.WSE.Vinit[iT]
    end

    return U, wind_vel
end

"""
    getWindSpeedT(::Velocity_Interpolation, wind_vel::AbstractMatrix, iT, t)

Interpolates the wind speed at a given time `t` using the provided wind velocity matrix `wind_vel`.
Uniform interpolation - all turbines experience the same changes.

# Arguments
- `::Velocity_Interpolation`: Specifies the interpolation strategy to use.
- `wind_vel::Matrix{Float64}`: A 2-column matrix: [time, wind_speed].
- `iT`: A single index or array of turbine indices.
- `t`: The requested time at which to interpolate the wind speed.

# Returns
- Interpolated wind speed at the respective turbine(s).
"""
function getWindSpeedT(::Velocity_Interpolation, wind_vel::AbstractMatrix, iT, t)
    times = wind_vel[:, 1]
    speeds = wind_vel[:, 2]

    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    itp = linear_interpolation(times, speeds, extrapolation_bc=Flat())
    u = itp(t)

    # Return a vector of wind speeds, one for each turbine index
    return fill(u, length(iT))
end

"""
    getWindSpeedT(::Velocity_Interpolation_wErrorCov, wind_vel::WindVelMatrix, iT, 
                  t)

Compute the wind speed at a given time `t` for the specified indices `iT` using the provided wind velocity data 
`wind_vel` and a velocity interpolation method with error covariance. 
Applies the same interpolated value across requested turbine indices.
Adds random noise from a given Cholesky decomposition matrix.

# Arguments
- `::Velocity_Interpolation_wErrorCov`: The interpolation method that includes error covariance handling.
- `wind_vel::WindVelMatrix`: See: [WindVelMatrix](@ref)
- `iT: Index or indices specifying which wind speed(s) to retrieve.
- `t`: The time at which to interpolate the wind speed.

# Returns
- The interpolated wind speed(s) at time `t` for the specified indices.
"""
function getWindSpeedT(::Velocity_Interpolation_wErrorCov, wind_vel::WindVelMatrix, iT, t)
    times = wind_vel.Data[:, 1]
    speeds = wind_vel.Data[:, 2]

    # Clamp time within data bounds
    if t < times[1]
        @warn "The time $t is out of bounds, using $(times[1]) instead."
        t_clamped = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, using $(times[end]) instead."
        t_clamped = times[end]
    else
        t_clamped = t
    end

    # Perform interpolation
    u = linear_interpolation(times, speeds, extrapolation_bc=Flat())(t_clamped)

    # Initialize velocities
    n = length(iT)
    Vel = fill(u, n)

    # Add randomness using Cholesky
    # Ensure noise dimensions match both n and CholSig dimensions
    chol_size = size(wind_vel.CholSig, 1)
    if n == chol_size
        noise = (wind_vel.CholSig * randn(RNG, chol_size))[1:n]
    else
        # If dimensions don't match, use appropriate subset or scaling
        if n <= chol_size
            noise = (wind_vel.CholSig * randn(RNG, chol_size))[1:n]
        else
            # If we need more noise values than CholSig size, repeat the pattern
            base_noise = wind_vel.CholSig * randn(RNG, chol_size)
            noise = repeat(base_noise, ceil(Int, n/chol_size))[1:n]
        end
    end
    Vel .= Vel .+ noise

    return Vel
end
"""
    getWindSpeedT(::Velocity_InterpTurbine, wind_vel::AbstractMatrix, iT, t)

Returns the wind speed at the specific turbine(s) and time.
The values are interpolated linearly between the set points.

# Arguments
- `::Velocity_InterpTurbine`: Marker for velocity interpolation at the turbine.
- `wind_vel::Matrix{Float64}`: Matrix where each row is time, U_T0, U_T1, ... U_Tn.
- `iT`: Index of the turbine for which the wind speed is requested.
- `t`: Time at which the wind speed is to be interpolated.

# Returns
- The interpolated wind speed at the specified turbine(s) and time as a `Float64`.
"""
function getWindSpeedT(::Velocity_InterpTurbine, wind_vel::AbstractMatrix, iT, t)
    times = wind_vel[:, 1]
    wind_data = wind_vel[:, 2:end]

    if t < times[1]So
        @warn "The time $t is out of bounds, using $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, using $(times[end]) instead."
        t = times[end]
    end

    # Interpolate specific turbine column(s)
    turbine_data = wind_data[:, iT]
    
    # Handle both single turbine (vector) and multiple turbines (matrix)
    if turbine_data isa AbstractVector
        # Single turbine: 1D interpolation
        itp = linear_interpolation(times, turbine_data, extrapolation_bc=Flat())
        return itp(t)
    else
        # Multiple turbines: 2D interpolation
        turbine_indices = collect(1:size(turbine_data, 2))
        itp = linear_interpolation((times, turbine_indices), turbine_data, extrapolation_bc=Flat())
        return [itp(t, i) for i in turbine_indices]
    end
end

"""
    getWindSpeedT(::Velocity_InterpTurbine_wErrorCov, wind_vel::WindVelMatrix, iT, t)

Returns the wind speed at the specified turbine(s) and time `t`,
including noise based on the covariance structure given in wind_vel.

- ::Velocity_InterpTurbine_wErrorCov: Model type indicating interpolation with covariance.
- wind_vel::WindVelMatrix: see [WindVelMatrix](@ref)
- iT: index or vector of indices of turbines
- t: time at which to query wind speed

Returns: Vector of wind speeds (with covariance noise) for iT at time t.
"""
function getWindSpeedT(::Velocity_InterpTurbine_wErrorCov, wind_vel::WindVelMatrix, iT, t)
    times = wind_vel.Data[:,1]
    speeds = wind_vel.Data[:,2:end]

    # Bounds check (clamp t)
    t_min = times[1]
    t_max = times[end]
    if t < t_min
        @warn "The time $t is out of bounds, will use $t_min instead."
        t = t_min
    elseif t > t_max
        @warn "The time $t is out of bounds, will use $t_max instead."
        t = t_max
    end

    # Create individual interpolants for each turbine column
    n_turbines = size(speeds, 2)
    interpolants = [interpolate((times,), speeds[:, i], Gridded(Linear())) for i in 1:n_turbines]

    # Evaluate all turbines at time `t`
    wind_vel_out = [itp(t) for itp in interpolants]

    vel = wind_vel_out[iT]               # select turbine(s) of interest

    # Add random error with given covariance (via Cholesky factor)
    # Simulate a noise vector (randn) multiplied by Cholesky factor
    noise = (wind_vel.CholSig * randn(RNG, length(vel)))
    vel_noisy = vel + noise
    return vel_noisy
end

# The following code cannot work, because the original Matlab function is not well defined.
#
# """
#     getWindSpeedT(::Velocity_RW_with_Mean, WindVelNow, WindVel::WindVelInitType)

# Calculates the wind speed at a given time step using a velocity model that includes both random walk and mean wind components.

# # Arguments
# - `::Velocity_RW_with_Mean`: The wind velocity model type indicating the use of random walk with mean wind.
# - `WindVelNow::WindVelInitType`: Initial wind velocity and covariance.
# - `WindVel`: The current wind velocity.

# # Returns
# - The computed wind speed at the current time step.
# """
# function getWindSpeedT(::Velocity_RW_with_Mean, WindVelNow::WindVelInitType, WindVel)
#     weightedRandN = randn(length(WindVelNow))
#     Vel = WindVelNow .+ (WindVel.CholSig' * weightedRandN) .+ WindVel.MeanPull .* (WindVel.Init .- WindVelNow)
# end

"""
    getWindSpeedT(::Velocity_ZOH_wErrorCov, vel::Vector, wind_vel_chol_sig::AbstractMatrix)

Computes the wind speed at a given time step using the zero-order hold (ZOH) method with error covariance.

# Arguments
- `::Velocity_ZOH_wErrorCov`: Type indicator for dispatch, representing the ZOH method with error covariance.
- `vel::Vector`: The velocity vector from the previous time step.
- `wind_vel_chol_sig::AbstractMatrix`: The Cholesky factor of the wind velocity error covariance matrix.

# Returns
- `wind_speed::Vector{Float64}`: The computed wind speed vector at the current time step.
"""
function getWindSpeedT(:: Velocity_ZOH_wErrorCov, vel::Vector, wind_vel_chol_sig::AbstractMatrix)
    # Generate standard normal random vector of same length as Vel
    noise = randn(RNG, length(vel))

    # Multiply by Cholesky factor to induce correlation
    correlated_noise = wind_vel_chol_sig * noise

    # Add the correlated noise to Vel
    vel .+= correlated_noise
end





