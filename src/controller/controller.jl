# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getYaw(::Yaw_SOWFA, ConYawData::AbstractMatrix, iT, t) -> Float64 or Vector{Float64}

Return the yaw angle at time `t` for the specified turbine(s) using linear interpolation.

# Arguments
- `::Yaw_SOWFA`: Controller type dispatch parameter for SOWFA-style yaw control
- `ConYawData::Matrix{Float64}`: Control data matrix where:
  - First column contains time values (in seconds)
  - Subsequent columns contain yaw angles for each turbine (in degrees)
- `iT`: Turbine index or indices to query:
  - `Integer`: Single turbine index (1-based)
  - `AbstractVector{<:Integer}`: Vector of turbine indices for multiple turbines
- `t::Real`: Requested time (in seconds)

# Returns
- `Float64`: Single yaw angle (in degrees) if `iT` is an integer
- `Vector{Float64}`: Vector of yaw angles (in degrees) if `iT` is a vector

# Behavior
- **Interpolation**: Uses linear interpolation between time points with flat extrapolation
- **Out-of-bounds handling**: If `t` is outside the time range, the function:
  - Issues a warning message
  - Clamps `t` to the nearest boundary (first or last time point)
- **Single time point**: If only one time point exists, returns the corresponding yaw value directly
- **Error handling**: Throws an error if `iT` is not an integer or vector of integers

# Data Format
The `con_yaw_data` matrix should have the structure:
```
[time₁  yaw₁₁  yaw₁₂  ...  yaw₁ₙ]
[time₂  yaw₂₁  yaw₂₂  ...  yaw₂ₙ]
[  ⋮      ⋮      ⋮    ⋱     ⋮  ]
[timeₘ  yawₘ₁  yawₘ₂  ...  yawₘₙ]
```
where `m` is the number of time steps and `n` is the number of turbines.

# Examples
```julia
# Example control data: 3 time points, 2 turbines
con_yaw_data = [0.0  10.0  5.0;   # t=0s: T1=10°, T2=5°
                1.0  15.0  10.0;  # t=1s: T1=15°, T2=10°
                2.0  20.0  15.0]  # t=2s: T1=20°, T2=15°

# Get yaw for turbine 1 at t=0.5s (interpolated)
yaw1 = getYaw(Yaw_SOWFA(), con_yaw_data, 1, 0.5)  # Returns 12.5°

# Get yaw for multiple turbines at t=1.5s
yaws = getYaw(Yaw_SOWFA(), con_yaw_data, [1, 2], 1.5)  # Returns [17.5°, 12.5°]

# Out-of-bounds time (will issue warning)
yaw_oob = getYaw(Yaw_SOWFA(), con_yaw_data, 1, 5.0)  # Returns 20.0° with warning
```

# See Also
- [`Yaw_SOWFA`](@ref): Controller type for SOWFA-style yaw control
- `Interpolations.linear_interpolation`: Underlying interpolation method used
"""
function getYaw(::Yaw_SOWFA, con_yaw_data::AbstractMatrix, iT, t)  
    if t < con_yaw_data[1, 1]
        @warn "The time $t is out of bounds, will use $(con_yaw_data[1, 1]) instead."
        t = con_yaw_data[1, 1]
    elseif t > con_yaw_data[end, 1]
        @warn "The time $t is out of bounds, will use $(con_yaw_data[end, 1]) instead."
        t = con_yaw_data[end, 1]
    end

    time = con_yaw_data[:, 1]
    yaw_data = con_yaw_data[:, 2:end]

    # Handle special case of single time point
    if length(time) == 1
        if isa(iT, Integer)
            return yaw_data[1, iT]
        elseif isa(iT, AbstractVector{<:Integer})
            return [yaw_data[1, i] for i in iT]
        else
            error("Invalid type for iT. Should be Integer or Vector of Integers.")
        end
    end

    # Create interpolation object for each turbine column
    interp_funcs = [linear_interpolation(time, yaw_data[:, j], extrapolation_bc=Flat()) for j in 1:size(yaw_data, 2)]

    # Get interpolated yaw(s)
    if isa(iT, Integer)
        return interp_funcs[iT](t)
    elseif isa(iT, AbstractVector{<:Integer})
        return [interp_funcs[i](t) for i in iT]
    else
        error("Invalid type for iT. Should be Integer or Vector of Integers.")
    end
end
