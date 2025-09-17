# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getYaw(::Yaw_SOWFA, con_yaw_data::AbstractMatrix, iT, t) -> Float64 or Vector{Float64}

Return the yaw angle at time `t` for the specified turbine(s) using linear interpolation.

# Arguments
- `::Yaw_SOWFA`: Controller type dispatch parameter for SOWFA-style yaw control
- `con_yaw_data::Matrix{Float64}`: Control data matrix where:
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

"""
    getYaw(::Yaw_Constant, con_yaw_data::AbstractMatrix, iT, t)

Return a single constant yaw angle (degrees) for one or multiple turbines.

Assumptions / Simplified Layout:
- `con_yaw_data` is a matrix with at least one row and one column.
- Only the value `con_yaw_data[1,1]` is used; any extra rows/columns are ignored.

Arguments:
- `::Yaw_Constant`: dispatch marker
- `con_yaw_data`: matrix whose first element holds the constant yaw angle
- `iT`: turbine index (Integer) or vector of indices
- `t`: ignored (kept for a uniform interface)

Returns:
- `Float64` when `iT` is an Integer
- `Vector{Float64}` when `iT` is a vector   
"""
function getYaw(::Yaw_Constant, con_yaw_data::AbstractMatrix, iT, t)
    size(con_yaw_data, 1) == 0 && error("con_yaw_data must have at least one row: $(con_yaw_data)")
    size(con_yaw_data, 2) == 0 && error("con_yaw_data must have at least one column: $(con_yaw_data)")
    yaw = con_yaw_data[1,1]
    if isa(iT, Integer)
        return yaw
    elseif isa(iT, AbstractVector{<:Integer})
        return fill(yaw, length(iT))
    else
        error("Invalid type for iT. Should be Integer or Vector of Integers.")
    end
end

"""
    getInduction(::Induction_Constant, con_induction_data::AbstractMatrix, iT, t) -> Float64 or Vector{Float64}

Return a single constant induction factor for one or multiple turbines.

# Arguments
- `::Induction_Constant`: Controller type dispatch parameter for constant induction control
- `con_induction_data::AbstractMatrix`: Control data matrix where only the first element is used
- `iT`: Turbine index or indices to query:
  - `Integer`: Single turbine index (1-based)
  - `AbstractVector{<:Integer}`: Vector of turbine indices for multiple turbines
- `t::Real`: Requested time (ignored for constant induction, kept for uniform interface)

# Returns
- `Float64`: Single induction factor if `iT` is an integer
- `Vector{Float64}`: Vector of induction factors if `iT` is a vector

# Behavior
- **Constant value**: Uses only `con_induction_data[1,1]` as the constant induction factor
- **Time independence**: The time parameter `t` is ignored since induction is constant
- **Multiple turbines**: Returns the same constant value for all requested turbines
- **Error handling**: Throws an error if `iT` is not an integer or vector of integers

# Data Format
The `con_induction_data` matrix should have at least one element:
- Only `con_induction_data[1,1]` is used
- Any additional rows/columns are ignored

# Examples
```julia
# Example control data with constant induction factor
con_induction_data = [0.33]  # Single constant induction factor

# Get induction for turbine 1 (any time)
induction1 = getInduction(Induction_Constant(), con_induction_data, 1, 100.0)  # Returns 0.33

# Get induction for multiple turbines
inductions = getInduction(Induction_Constant(), con_induction_data, [1, 2, 3], 50.0)  # Returns [0.33, 0.33, 0.33]
```

# See Also
- [`Induction_Constant`](@ref): Controller type for constant induction control
"""
function getInduction(::Induction_Constant, con::Con, iT, t)
    induction = con.induction_fixed
    if isa(iT, Integer)
        return induction
    elseif isa(iT, AbstractVector{<:Integer})
        return fill(induction, length(iT))
    else
        error("Invalid type for iT. Should be Integer or Vector of Integers.")
    end
end

"""
    getInduction(::Induction_MPC, con_induction_data::AbstractMatrix, iT, t) -> Float64 or Vector{Float64}

Return the induction factor at time `t` for the specified turbine(s) using linear interpolation.

# Arguments
- `::Induction_MPC`: Controller type dispatch parameter for MPC-style induction control
- `con_induction_data::Matrix{Float64}`: Control data matrix where:
  - First column contains time values (in seconds)
  - Subsequent columns contain induction factors for each turbine (dimensionless)
- `iT`: Turbine index or indices to query:
  - `Integer`: Single turbine index (1-based)
  - `AbstractVector{<:Integer}`: Vector of turbine indices for multiple turbines
- `t::Real`: Requested time (in seconds)

# Returns
- `Float64`: Single induction factor if `iT` is an integer
- `Vector{Float64}`: Vector of induction factors if `iT` is a vector

# Behavior
- **Interpolation**: Uses linear interpolation between time points with flat extrapolation
- **Out-of-bounds handling**: If `t` is outside the time range, the function:
  - Issues a warning message
  - Clamps `t` to the nearest boundary (first or last time point)
- **Single time point**: If only one time point exists, returns the corresponding induction value directly
- **Error handling**: Throws an error if `iT` is not an integer or vector of integers

# Data Format
The `con_induction_data` matrix should have the structure:
```
[time₁  induction₁₁  induction₁₂  ...  induction₁ₙ]
[time₂  induction₂₁  induction₂₂  ...  induction₂ₙ]
[  ⋮        ⋮            ⋮        ⋱       ⋮     ]
[timeₘ  inductionₘ₁  inductionₘ₂  ...  inductionₘₙ]
```
where `m` is the number of time steps and `n` is the number of turbines.

# Examples
```julia
# Example control data: 3 time points, 2 turbines
con_induction_data = [0.0  0.30  0.25;   # t=0s: T1=0.30, T2=0.25
                      1.0  0.35  0.30;   # t=1s: T1=0.35, T2=0.30
                      2.0  0.40  0.35]   # t=2s: T1=0.40, T2=0.35

# Get induction for turbine 1 at t=0.5s (interpolated)
induction1 = getInduction(Induction_MPC(), con_induction_data, 1, 0.5)  # Returns 0.325

# Get induction for multiple turbines at t=1.5s
inductions = getInduction(Induction_MPC(), con_induction_data, [1, 2], 1.5)  # Returns [0.375, 0.325]

# Out-of-bounds time (will issue warning)
induction_oob = getInduction(Induction_MPC(), con_induction_data, 1, 5.0)  # Returns 0.40 with warning
```

# See Also
- [`Induction_MPC`](@ref): Controller type for MPC-style induction control
- `Interpolations.linear_interpolation`: Underlying interpolation method used
"""
function getInduction(::Induction_MPC, con::Con, iT, t) 
    con_induction_data = con.induction_data
    if t < con_induction_data[1, 1]
        @warn "The time $t is out of bounds, will use $(con_induction_data[1, 1]) instead."
        t = con_induction_data[1, 1]
    elseif t > con_induction_data[end, 1]
        @warn "The time $t is out of bounds, will use $(con_induction_data[end, 1]) instead."
        t = con_induction_data[end, 1]
    end

    time = con_induction_data[:, 1]
    induction_data = con_induction_data[:, 2:end]

    # Handle special case of single time point
    if length(time) == 1
        if isa(iT, Integer)
            return induction_data[1, iT]
        elseif isa(iT, AbstractVector{<:Integer})
            return [induction_data[1, i] for i in iT]
        else
            error("Invalid type for iT. Should be Integer or Vector of Integers.")
        end
    end

    # Create interpolation object for each turbine column
    interp_funcs = [linear_interpolation(time, induction_data[:, j], extrapolation_bc=Flat()) for j in 1:size(induction_data, 2)]

    # Get interpolated induction(s)
    if isa(iT, Integer)
        return interp_funcs[iT](t)
    elseif isa(iT, AbstractVector{<:Integer})
        return [interp_funcs[i](t) for i in iT]
    else
        error("Invalid type for iT. Should be Integer or Vector of Integers.")
    end
end
