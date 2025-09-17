# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getYaw(::Yaw_SOWFA, con::Con, iT, t) -> Float64 or Vector{Float64}

Return the yaw angle at time `t` for the specified turbine(s) using linear interpolation from controller data.

# Arguments
- `::Yaw_SOWFA`: Controller type dispatch parameter for SOWFA-style yaw control
- `con::Con`: Controller configuration struct containing yaw control data
- `iT`: Turbine index or indices to query:
  - `Integer`: Single turbine index (1-based)
  - `AbstractVector{<:Integer}`: Vector of turbine indices for multiple turbines
- `t::Real`: Requested time (in seconds)

# Returns
- `Float64`: Single yaw angle (in degrees) if `iT` is an integer
- `Vector{Float64}`: Vector of yaw angles (in degrees) if `iT` is a vector

# Behavior
- **Interpolation**: Uses linear interpolation between time points with flat extrapolation
- **Data source**: Retrieves yaw data from `con.yaw_data` matrix
- **Out-of-bounds handling**: If `t` is outside the time range, the function:
  - Issues a warning message
  - Clamps `t` to the nearest boundary (first or last time point)
- **Single time point**: If only one time point exists, returns the corresponding yaw value directly
- **Error handling**: Throws an error if `iT` is not an integer or vector of integers

# Data Format
The `con.yaw_data` matrix should have the structure:
```
[time₁  yaw₁₁  yaw₁₂  ...  yaw₁ₙ]
[time₂  yaw₂₁  yaw₂₂  ...  yaw₂ₙ]
[  ⋮      ⋮      ⋮    ⋱     ⋮  ]
[timeₘ  yawₘ₁  yawₘ₂  ...  yawₘₙ]
```
where `m` is the number of time steps and `n` is the number of turbines.

# Examples
```julia
# Create controller configuration with yaw data
con = Con(yaw="SOWFA", yaw_data=[0.0  10.0  5.0;   # t=0s: T1=10°, T2=5°
                                 1.0  15.0  10.0;  # t=1s: T1=15°, T2=10°
                                 2.0  20.0  15.0]) # t=2s: T1=20°, T2=15°

# Get yaw for turbine 1 at t=0.5s (interpolated)
yaw1 = getYaw(Yaw_SOWFA(), con, 1, 0.5)  # Returns 12.5°

# Get yaw for multiple turbines at t=1.5s
yaws = getYaw(Yaw_SOWFA(), con, [1, 2], 1.5)  # Returns [17.5°, 12.5°]

# Out-of-bounds time (will issue warning)
yaw_oob = getYaw(Yaw_SOWFA(), con, 1, 5.0)  # Returns 20.0° with warning
```

# See Also
- [`Yaw_SOWFA`](@ref): Controller type for SOWFA-style yaw control
- [`Con`](@ref): Controller configuration struct
- `Interpolations.linear_interpolation`: Underlying interpolation method used
"""
function getYaw(::Yaw_SOWFA, con::Con, iT, t)  
    con_yaw_data = con.yaw_data
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
    getYaw(::Yaw_Constant, con::Con, iT, t) -> Float64 or Vector{Float64}

Return a single constant yaw angle for one or multiple turbines.

# Arguments
- `::Yaw_Constant`: Controller type dispatch parameter for constant yaw control
- `con::Con`: Controller configuration struct containing fixed yaw value
- `iT`: Turbine index or indices to query:
  - `Integer`: Single turbine index (1-based) 
  - `AbstractVector{<:Integer}`: Vector of turbine indices for multiple turbines
- `t::Real`: Requested time (ignored for constant yaw, kept for uniform interface)

# Returns
- `Float64`: Single yaw angle (in degrees) if `iT` is an integer
- `Vector{Float64}`: Vector of yaw angles (in degrees) if `iT` is a vector

# Behavior
- **Constant value**: Uses `con.yaw_fixed` as the yaw angle for all turbines and times
- **Time independence**: The time parameter `t` is ignored since yaw is constant
- **Multiple turbines**: Returns the same constant value for all requested turbines
- **Error handling**: Throws an error if `iT` is not an integer or vector of integers

# Examples
```julia
# Create controller configuration with constant yaw
con = Con(yaw="Constant", yaw_fixed=270.0)  # 270° constant yaw

# Get yaw for turbine 1 (any time)
yaw1 = getYaw(Yaw_Constant(), con, 1, 100.0)  # Returns 270.0°

# Get yaw for multiple turbines
yaws = getYaw(Yaw_Constant(), con, [1, 2, 3], 50.0)  # Returns [270.0°, 270.0°, 270.0°]
```

# See Also
- [`Yaw_Constant`](@ref): Controller type for constant yaw control
- [`Con`](@ref): Controller configuration struct
"""
function getYaw(::Yaw_Constant, con::Con, iT, t)
    yaw = con.yaw_fixed
    if isa(iT, Integer)
        return yaw
    elseif isa(iT, AbstractVector{<:Integer})
        return fill(yaw, length(iT))
    else
        error("Invalid type for iT. Should be Integer or Vector of Integers.")
    end
end

"""
    getInduction(::Induction_Constant, con::Con, iT, t) -> Float64 or Vector{Float64}

Return a single constant induction factor for one or multiple turbines.

# Arguments
- `::Induction_Constant`: Controller type dispatch parameter for constant induction control
- `con::Con`: Controller configuration struct containing fixed induction value
- `iT`: Turbine index or indices to query:
  - `Integer`: Single turbine index (1-based)
  - `AbstractVector{<:Integer}`: Vector of turbine indices for multiple turbines
- `t::Real`: Requested time (ignored for constant induction, kept for uniform interface)

# Returns
- `Float64`: Single induction factor if `iT` is an integer
- `Vector{Float64}`: Vector of induction factors if `iT` is a vector

# Behavior
- **Constant value**: Uses `con.induction_fixed` as the induction factor for all turbines and times
- **Time independence**: The time parameter `t` is ignored since induction is constant
- **Multiple turbines**: Returns the same constant value for all requested turbines
- **Error handling**: Throws an error if `iT` is not an integer or vector of integers

# Examples
```julia
# Create controller configuration with constant induction
con = Con(induction="Constant", induction_fixed=0.33)  # 0.33 constant induction factor

# Get induction for turbine 1 (any time)
induction1 = getInduction(Induction_Constant(), con, 1, 100.0)  # Returns 0.33

# Get induction for multiple turbines
inductions = getInduction(Induction_Constant(), con, [1, 2, 3], 50.0)  # Returns [0.33, 0.33, 0.33]
```

# See Also
- [`Induction_Constant`](@ref): Controller type for constant induction control
- [`Con`](@ref): Controller configuration struct
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
    getInduction(::Induction_MPC, con::Con, iT, t) -> Float64 or Vector{Float64}

Return the induction factor at time `t` for the specified turbine(s) using linear interpolation from controller data.

# Arguments
- `::Induction_MPC`: Controller type dispatch parameter for MPC-style induction control
- `con::Con`: Controller configuration struct containing induction control data
- `iT`: Turbine index or indices to query:
  - `Integer`: Single turbine index (1-based)
  - `AbstractVector{<:Integer}`: Vector of turbine indices for multiple turbines
- `t::Real`: Requested time (in seconds)

# Returns
- `Float64`: Single induction factor if `iT` is an integer
- `Vector{Float64}`: Vector of induction factors if `iT` is a vector

# Behavior
- **Interpolation**: Uses linear interpolation between time points with flat extrapolation
- **Data source**: Retrieves induction data from `con.induction_data` matrix
- **Out-of-bounds handling**: If `t` is outside the time range, the function:
  - Issues a warning message
  - Clamps `t` to the nearest boundary (first or last time point)
- **Single time point**: If only one time point exists, returns the corresponding induction value directly
- **Error handling**: Throws an error if `iT` is not an integer or vector of integers

# Data Format
The `con.induction_data` matrix should have the structure:
```
[time₁  induction₁₁  induction₁₂  ...  induction₁ₙ]
[time₂  induction₂₁  induction₂₂  ...  induction₂ₙ]
[  ⋮        ⋮            ⋮        ⋱       ⋮     ]
[timeₘ  inductionₘ₁  inductionₘ₂  ...  inductionₘₙ]
```
where `m` is the number of time steps and `n` is the number of turbines.

# Examples
```julia
# Create controller configuration with induction data
con = Con(induction="MPC", 
          induction_data=[0.0  0.30  0.25;   # t=0s: T1=0.30, T2=0.25
                         1.0  0.35  0.30;   # t=1s: T1=0.35, T2=0.30
                         2.0  0.40  0.35])  # t=2s: T1=0.40, T2=0.35

# Get induction for turbine 1 at t=0.5s (interpolated)
induction1 = getInduction(Induction_MPC(), con, 1, 0.5)  # Returns 0.325

# Get induction for multiple turbines at t=1.5s
inductions = getInduction(Induction_MPC(), con, [1, 2], 1.5)  # Returns [0.375, 0.325]

# Out-of-bounds time (will issue warning)
induction_oob = getInduction(Induction_MPC(), con, 1, 5.0)  # Returns 0.40 with warning
```

# See Also
- [`Induction_MPC`](@ref): Controller type for MPC-style induction control
- [`Con`](@ref): Controller configuration struct
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
