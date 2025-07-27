# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getYaw(::Yaw_SOWFA, ConYawData::Matrix{Float64}, iT, t)

Return the yaw angle(s) at time `t` for specified turbine index/indices `iT` using SOWFA yaw controller data.

# Arguments
- `::Yaw_SOWFA`: Controller model type indicating SOWFA-based yaw control
- `ConYawData::Matrix{Float64}`: Matrix where the first column contains time values and subsequent 
  columns contain yaw angles for each turbine
- `iT`: Turbine index (Integer) or vector of turbine indices (AbstractVector{<:Integer})
- `t`: Requested time value for yaw angle interpolation

# Returns
- For single turbine index: Float64 representing the yaw angle at time `t`
- For multiple turbine indices: Vector{Float64} containing yaw angles for all requested turbines at time `t`

# Description
This function performs temporal interpolation of yaw angle data from SOWFA simulations. It uses linear 
interpolation between time points and applies flat extrapolation for times outside the data range.
The function automatically handles boundary cases and issues warnings when the requested time is 
outside the available data range.

# Examples
```julia
# Single turbine yaw at time t=100.0
yaw_angle = getYaw(Yaw_SOWFA(), yaw_data_matrix, 1, 100.0)

# Multiple turbines yaw at time t=100.0
yaw_angles = getYaw(Yaw_SOWFA(), yaw_data_matrix, [1, 2, 3], 100.0)
```

# Notes
- Time values outside the data range will trigger a warning and use boundary values
- For single time point datasets, direct lookup is performed without interpolation
- Uses linear interpolation with flat extrapolation for robust temporal estimation

# See also
- [`Yaw_SOWFA`](@ref): SOWFA-based yaw controller model type
"""
function getYaw(::Yaw_SOWFA, ConYawData::Matrix{Float64}, iT, t)
    if t < ConYawData[1, 1]
        @warn "The time $t is out of bounds, will use $(ConYawData[1, 1]) instead."
        t = ConYawData[1, 1]
    elseif t > ConYawData[end, 1]
        @warn "The time $t is out of bounds, will use $(ConYawData[end, 1]) instead."
        t = ConYawData[end, 1]
    end

    time = ConYawData[:, 1]
    yaw_data = ConYawData[:, 2:end]

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
