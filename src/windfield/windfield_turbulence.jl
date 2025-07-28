# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getWindTiT(::TI_Constant, wind_ti::Number, iT, _)

Return turbulence intensity for the requested turbine(s).

# Arguments
- `::TI_Constant`: type parameter to indicate constant wind turbulence
- `wind_ti::Number`: Constant value (turbulence intensity)
- `iT`: Index or indices of the turbines
- `_`: will be ignored

# Returns
- `Ti`: Array of turbulence intensity values for each turbine index
"""
function getWindTiT(::TI_Constant, wind_ti, iT, _)
    if isa(iT, AbstractArray)
        return fill(wind_ti, size(iT))
    else
        return wind_ti
    end
end

# function Ti = getWindTiT(wind_ti,iT,t)
# %GETWINDTI Return turbulence intensity for the requested turbine(s)
# % Uniform interpolation version - all turbines experience the same changes
# %   Expects a .csv called "WindTI" with the structure
# %       time TI
# %       time TI
# %        ...
# %       time TI
# %   The value is interpolated linearly between the setpoints
# % ======= Input ======
# % wind_ti    = (t,TI) pairs between which is linearly interpolated
# % iT        = Index/Indeces of the turbines
# % t         = time of request
# % ======================================================================= %

# if t<wind_ti(1,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(wind_ti(1,1)) ' instead.']);
#     t = wind_ti(1,1);
# elseif t>wind_ti(end,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(wind_ti(end,1)) ' instead.']);
#     t = wind_ti(end,1);
# end

# Ti = ones(size(iT))*interp1(wind_ti(:,1),wind_ti(:,2),t);
# end

"""
    getWindTiT(::TI_Interpolation, wind_ti::AbstractMatrix, iT, t)

Interpolates the wind turbulence intensity (TI) at a given time `t` using the specified `TI_Interpolation` method.

# Arguments
- `::TI_Interpolation`: Use linear interpolation to calculate the turbulence intensity.
- `wind_ti::Matrix`: Matrix containing wind turbulence intensity values over time.
- `iT`: Index/indices of the turbines (can be Int or array).
- `t`: The specific time at which to interpolate the turbulence intensity.

# Returns
- The interpolated turbulence for the requested turbine(s) at time `t`.

# Notes
- The function assumes that `wind_ti` contains the necessary data for interpolation
  as (time, TI) pairs (n×2 matrix)
- Uniform interpolation version - all turbines experience the same changes.
"""
function getWindTiT(::TI_Interpolation, wind_ti::AbstractMatrix, iT, t)

    # Extract time and TI columns
    times = wind_ti[:, 1]
    TIs = wind_ti[:, 2]

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
    if isa(iT, AbstractArray)
        return fill(Ti_val, size(iT))
    else
        return Ti_val
    end
end

# function Ti = getWindTiT(wind_ti,iT,t)
# %GETWINDTI Return turbulence intensity for the requested turbine(s)
# % ======================================================================= %
# % Individual Turbine value implementation
# %   requires a .csv in the simulation folder called WindTITurbine.csv
# %   where each row is a
# %       time, TI_T0, TI_T1, ... TI_Tn
# %   setpoint in time. The values are interploated linearly between the
# %   setpoints.
# % ======= Input ======
# % wind_ti    = (t,TI_T0, TI_T1, ... TI_Tn)
# % iT        = Index/Indeces of the turbines
# % t         = time of request
# % ======================================================================= %

# if t<wind_ti(1,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(wind_ti(1,1)) ' instead.']);
#     t = wind_ti(1,1);
# elseif t>wind_ti(end,1)
#     warning(['The time ' num2str(t) ' is out of bounds, will use '...
#         num2str(wind_ti(end,1)) ' instead.']);
#     t = wind_ti(end,1);
# end
# Ti_out = interp1(wind_ti(:,1),wind_ti(:,2:end),t);
# Ti = Ti_out(iT);
# end

"""
    getWindTiT(::TI_InterpTurbine, wind_ti::AbstractMatrix, iT, t)

Retrieve the wind turbulence intensity (TI) for a specific turbine at a given time.

# Arguments
- `::TI_InterpTurbine`: The turbulence intensity interpolation object for the turbine.
- `wind_ti::AbstractMatrix`: Matrix containing wind turbulence intensity values.
- `iT`: Index of the turbine for which the TI is requested.
- `t`: Time at which the TI value is needed.

# Returns
- The interpolated wind turbulence intensity value for the specified turbine at time `t`.
"""
function getWindTiT(::TI_InterpTurbine, wind_ti::AbstractMatrix, iT, t)
    # wind_ti: Matrix (time, TI_T0, TI_T1, ..., TI_Tn)
    # iT: Index or indices of turbines (1-based)
    # t: Requested time

    times = wind_ti[:, 1]
    n_turbines = size(wind_ti, 2) - 1

    # Clamp t to the bounds of times
    if t < times[1]
        @warn "The time $t is out of bounds, will use $(times[1]) instead."
        t = times[1]
    elseif t > times[end]
        @warn "The time $t is out of bounds, will use $(times[end]) instead."
        t = times[end]
    end

    # Interpolate for each turbine
    Ti_out = [interp1d(times, wind_ti[:, j+1], t) for j in 1:n_turbines]

    # Return value(s) for the requested turbine(s)
    return Ti_out[iT]
end

# Helper function for linear interpolation at a single point
function interp1d(x, y, xi)
    # Handle the case where all y values are the same (constant)
    if all(y .== y[1])
        return y[1]
    end
    
    # Handle case with only one data point
    if length(x) == 1
        return y[1]
    end
    
    itp = linear_interpolation(x, y, extrapolation_bc=Flat())
    return itp(xi)
end


