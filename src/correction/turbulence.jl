# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    correctTI!(::TI_None, set, wf, wind, sim_time)

Update turbulence intensity values in the wind farm state matrix without correction.

This function implements the "no correction" strategy for turbulence intensity, where 
the wind farm turbulence intensity values are updated with fresh data from the wind 
field model without applying any correction algorithms. It serves as the baseline 
approach for turbulence intensity handling in FLORIDyn simulations.

# Arguments
- `::TI_None`: Dispatch type indicating no turbulence intensity correction algorithm
- `set`: Settings object containing simulation configuration and turbulence model parameters
- `wf`: Wind farm object containing the state matrices to be updated
  - `wf.States_WF`: Wind field states matrix where column 3 contains turbulence intensity values
  - `wf.StartI`: Starting indices for each turbine's operational points
  - `wf.nT`: Number of turbines
- `wind`: Wind configuration object containing turbulence intensity input specifications
  - `wind.input_ti`: String specifying the turbulence intensity input method
  - `wind.ti`: Turbulence intensity data or model parameters
- `sim_time`: Current simulation time for time-dependent turbulence intensity retrieval

# Returns
- `nothing`: The function modifies the wind farm state in-place

# Notes
- The function modifies the wind farm object in-place (indicated by the `!` suffix)
- This "no correction" approach provides baseline turbulence intensity without 
  applying wake-induced corrections or measurement-based adjustments
"""
function correctTI!(::TI_None, set, wf, wind, sim_time)
    # correctTI! updates the turbulent intensity (TI) value inT.States_WF
    # at the rowT.StartI and column 3, using data from getDataTI.

    # Get new turbulent intensity value
    TI = getDataTI(set, wind, wf, sim_time)

    # Update the TI value in the States_WF matrix
    try
       wf.States_WF[wf.StartI, 3] .= TI'
    catch e
        error("Error updating wf.States_WF: $(e.msg)")
    end

    return nothing
end

"""
    getDataTI(set, wind, wf, sim_time)

Retrieve turbulence intensity data based on the specified input method.

This function serves as the dispatcher for obtaining turbulence intensity 
values from various sources and models. It supports multiple input methods including 
Ensemble Kalman Filter approaches, zero-order hold strategies, random walk models, 
and closed-loop control with weighted interpolation.

# Arguments
- `set`: Settings object containing simulation configuration
  - `set.turb_mode`: Turbulence model configuration for standard retrieval methods
- `wind`: Wind configuration object specifying the turbulence intensity source
  - `wind.input_ti`: String specifying the input method for turbulence intensity
  - `wind.ti`: Turbulence intensity data, model parameters, or covariance matrices
- `wf`: Wind farm object containing current simulation state
  - `wf.States_WF`: Wind field states matrix with turbulence intensity in column 3
  - `wf.StartI`: Starting indices for each turbine's operational points  
  - `wf.nT`: Number of turbines
  - `wf.C_TI`: Interpolation coefficient matrix (for weighted methods)
- `sim_time`: Current simulation time for time-dependent turbulence intensity models

# Returns
- `TI`: Turbulence intensity values for all turbines as a vector or array

# Notes
- Random walk perturbations use Cholesky decomposition for proper covariance
- Interpolation methods require precomputed coefficient matrices
- Time-dependent methods support dynamic turbulence intensity evolution
"""
function getDataTI(set, wind, wf, sim_time)
    # GETDATATI retrieves the data for the ambient turbulence intensity

    if wind.input_ti == "EnKF_InterpTurbine"
        TI = getWindTiT_EnKF(set.turb_mode, wind.ti, collect(1:wf.nT), sim_time)
    elseif wind.input_ti == "EnKF_ZOH"
        TI =wf.States_WF[wf.StartI, 3]
    elseif wind.input_ti == "EnKF_RW"
        TI =wf.States_WF[wf.StartI, 3]
        # 'phi' is not defined in your snippet. Replace with wf.nT if that's correct.
        TI = TI .+ (randn(1, wf.nT) * wind.ti.CholSig)'
    elseif wind.input_ti == "CLC_weighted_ZOH"
        TI = wf.C_TI *wf.States_WF[:, 3]
    else
        TI = getWindTiT(set.turb_mode, wind.ti, collect(1:wf.nT), sim_time)
    end

    return TI
end

