# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function getYaw(::Yaw_SOWFA, ConYawData::Matrix{Float64}, iT, t)
    # GETYAW Return the yaw angle at time t for turbine index/indices iT
    # ConYawData: matrix where first column is time, and subsequent columns are yaw angles for turbines
    # iT: single index or array of turbine indices (1-based indexing)
    # t: requested time
    
    if t < ConYawData[1, 1]
        @warn "The time $t is out of bounds, will use $(ConYawData[1, 1]) instead."
        t = ConYawData[1, 1]
    elseif t > ConYawData[end, 1]
        @warn "The time $t is out of bounds, will use $(ConYawData[end, 1]) instead."
        t = ConYawData[end, 1]
    end

    time = ConYawData[:, 1]
    yaw_data = ConYawData[:, 2:end]

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
