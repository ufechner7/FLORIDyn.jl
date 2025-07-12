# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getWindTiT(::TI_Constant, WindTi, iT)

Return turbulence intensity for the requested turbine(s).

# Arguments
- `::TI_Constant`: type parameter to indicate constant wind turbulence
- `WindTi`: Constant value (turbulence intensity)
- `iT`: Index or indices of the turbines

# Returns
- `Ti`: Array of turbulence intensity values for each turbine index
"""
function getWindTiT(::TI_Constant, WindTi, iT)
    return fill(WindTi, size(iT))
end

# function Ti = getWindTiT_EnKF(WindTi,iT,t)
# %GETWINDTI Return turbulence intensity for the requested turbine(s)
# % ======================================================================= %
# % Individual Turbine value implementation
# %   requires a .csv in the simulation folder called WindTITurbine.csv
# %   where each row is a
# %       time, TI_T0, TI_T1, ... TI_Tn
# %   setpoint in time. The values are interploated linearly between the
# %   setpoints.
# % ======= Input ======
# % WindTi    = (t,TI_T0, TI_T1, ... TI_Tn)
# % iT        = Index/Indeces of the turbines
# % t         = time of request
# % ======================================================================= %

# if t<WindTi(1,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(WindTi(1,1)) ' instead.']);
#     t = WindTi(1,1);
# elseif t>WindTi(end,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(WindTi(end,1)) ' instead.']);
#     t = WindTi(end,1);
# end
# Ti_out = interp1(WindTi(:,1),WindTi(:,2:end),t);
# Ti = Ti_out(iT);
# end

"""
    getWindTiT_EnKF(::TI_EnKF_InterpTurbine, WindTi, iT, t)

Return turbulence intensity for the requested turbine(s) at time `t`.

# Arguments
- `::TI_EnKF_InterpTurbine`: type parameter to indicate Individual Turbine value implementation
- `WindTi`: Matrix where each row is [time, TI_T0, TI_T1, ... TI_Tn]
- `iT`: Index or indices of the turbines (1-based)
- `t`: Time of request

# Returns
- `Ti`: Turbulence intensity for the requested turbines at time `t`
"""
function getWindTiT_EnKF(::TI_EnKF_InterpTurbine, WindTi, iT, t)
    times = WindTi[:, 1]
    TI = WindTi[:, 2:end]

    # Clamp t to the bounds of available times
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Interpolate each turbine column at time t
    Ti_out = [linear_interpolation(times, TI[:, j], extrapolation_bc=Line())(t) for j in 1:size(TI, 2)]
    Ti_out[iT]
end

# function Ti = getWindTiT(WindTi,iT,t)
# %GETWINDTI Return turbulence intensity for the requested turbine(s)
# % Uniform interpolation version - all turbines experience the same changes
# %   Expects a .csv called "WindTI" with the structure
# %       time TI
# %       time TI
# %        ...
# %       time TI
# %   The value is interpolated linearly between the setpoints
# % ======= Input ======
# % WindTi    = (t,TI) pairs between which is linearly interpolated
# % iT        = Index/Indeces of the turbines
# % t         = time of request
# % ======================================================================= %

# if t<WindTi(1,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(WindTi(1,1)) ' instead.']);
#     t = WindTi(1,1);
# elseif t>WindTi(end,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(WindTi(end,1)) ' instead.']);
#     t = WindTi(end,1);
# end

# Ti = ones(size(iT))*interp1(WindTi(:,1),WindTi(:,2),t);
# end

"""
    getWindTiT(::TI_Interpolation, WindTi::Matrix{<:Real}, iT, t::Real)

Interpolates the wind turbulence intensity (TI) at a given time `t` using the specified `TI_Interpolation` method.

# Arguments
- `::TI_Interpolation`: Use linear interpolation to calculate the turbulence intensity.
- `WindTi::Matrix{<:Real}`: Matrix containing wind turbulence intensity values over time.
- `iT`: Index/indices of the turbines (can be Int or array).
- `t::Real`: The specific time at which to interpolate the turbulence intensity.

# Returns
- The interpolated turbulence for the requested turbine(s) at time `t`.

# Notes
- The function assumes that `WindTi` contains the necessary data for interpolation
  as (time, TI) pairs (n×2 matrix)
- Uniform interpolation version - all turbines experience the same changes.
"""
function getWindTiT(::TI_Interpolation, WindTi::Matrix{<:Real}, iT, t::Real)

    # Extract time and TI columns
    times = WindTi[:, 1]
    TIs = WindTi[:, 2]

    # Clamp t to the bounds, with warnings
    if t < times[1]
        @warn("The time $t is out of bounds, will use $(times[1]) instead.")
        t = times[1]
    elseif t > times[end]
        @warn("The time $t is out of bounds, will use $(times[end]) instead.")
        t = times[end]
    end

    # Linear interpolation
    # Search for the interval
    idx = searchsortedlast(times, t)
    if idx == length(times)
        Ti_val = TIs[end]
    elseif times[idx] == t
        Ti_val = TIs[idx]
    else
        # Linear interpolation
        t0, t1 = times[idx], times[idx+1]
        ti0, ti1 = TIs[idx], TIs[idx+1]
        Ti_val = ti0 + (ti1 - ti0) * (t - t0) / (t1 - t0)
    end

    # Output: ones(size(iT)) * Ti_val
    # If iT is a single Int, return scalar; if array, return array
    fill(Ti_val, length(iT))
end

# function Ti = getWindTiT(WindTi,iT,t)
# %GETWINDTI Return turbulence intensity for the requested turbine(s)
# % ======================================================================= %
# % Individual Turbine value implementation
# %   requires a .csv in the simulation folder called WindTITurbine.csv
# %   where each row is a
# %       time, TI_T0, TI_T1, ... TI_Tn
# %   setpoint in time. The values are interploated linearly between the
# %   setpoints.
# % ======= Input ======
# % WindTi    = (t,TI_T0, TI_T1, ... TI_Tn)
# % iT        = Index/Indeces of the turbines
# % t         = time of request
# % ======================================================================= %

# if t<WindTi(1,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(WindTi(1,1)) ' instead.']);
#     t = WindTi(1,1);
# elseif t>WindTi(end,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(WindTi(end,1)) ' instead.']);
#     t = WindTi(end,1);
# end
# Ti_out = interp1(WindTi(:,1),WindTi(:,2:end),t);
# Ti = Ti_out(iT);
# end

"""
    getWindTiT(::TI_InterpTurbine, WindTi, iT, t)

Retrieve the wind turbulence intensity (TI) for a specific turbine at a given time.

# Arguments
- `::TI_InterpTurbine`: The turbulence intensity interpolation object for the turbine.
- `WindTi`: Array or data structure containing wind turbulence intensity values.
- `iT`: Index of the turbine for which the TI is requested.
- `t`: Time at which the TI value is needed.

# Returns
- The interpolated wind turbulence intensity value for the specified turbine at time `t`.
"""
function getWindTiT(::TI_InterpTurbine, WindTi, iT, t)
    # WindTi: Matrix (time, TI_T0, TI_T1, ..., TI_Tn)
    # iT: Index or indices of turbines (1-based)
    # t: Requested time

    times = WindTi[:, 1]
    n_turbines = size(WindTi, 2) - 1

    # Clamp t to the bounds of times
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Interpolate for each turbine
    Ti_out = [interp1d(times, WindTi[:, j+1], t) for j in 1:n_turbines]

    # Return value(s) for the requested turbine(s)
    return Ti_out[iT]
end

# Helper function for linear interpolation at a single point
function interp1d(x, y, xi)
    itp = LinearInterpolation(x, y, extrapolation_bc=Flat())
    return itp(xi)
end


