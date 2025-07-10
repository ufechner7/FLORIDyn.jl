# function Ti = getWindTiT(WindTi,iT,~)
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
