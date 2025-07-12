# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    getWindTiT(::TI_Constant, WindTi, iT)

Return turbulence intensity for the requested turbine(s).

# Arguments
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
