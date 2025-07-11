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
