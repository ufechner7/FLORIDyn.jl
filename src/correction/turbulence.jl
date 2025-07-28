# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    correctTI!(::TI_None, set::Settings, wf::WindFarm, wind::Wind, t) -> Nothing

Update turbulence intensity values in the wind farm state matrix without correction.

This function implements the "no correction" strategy for turbulence intensity, where 
the wind farm turbulence intensity values are updated with fresh data from the wind 
field model without applying any correction algorithms. It serves as the baseline 
approach for turbulence intensity handling in FLORIDyn simulations.

# Arguments
- `::TI_None`: Dispatch type indicating no turbulence intensity correction algorithm
- `set::Settings`: Settings object containing simulation configuration and turbulence model parameters
  - `set.turb_mode`: Turbulence model configuration specifying the retrieval method
- `wf::WindFarm`: Wind farm object containing the state matrices to be updated
  - `wf.States_WF`: Wind field states matrix where column 3 contains turbulence intensity values
  - `wf.StartI`: Starting indices for each turbine's operational points
  - `wf.nT`: Number of turbines
- `wind::Wind`: Wind configuration object containing turbulence intensity data
  - `wind.ti`: Turbulence intensity data or model parameters
- `t`: Current simulation time for time-dependent turbulence intensity retrieval

# Returns
- `Nothing`: The function modifies the wind farm state in-place

# Behavior
1. Retrieves current turbulence intensity values for all turbines using [`getDataTI`](@ref)
2. Updates the wind farm state matrix `wf.States_WF` at rows `wf.StartI` and column 3
3. Transposes the turbulence intensity vector to match the matrix structure
4. Provides error handling for matrix update operations

# Example
```julia
# Update turbulence intensity without correction at t=50s
correctTI!(TI_None(), settings, wf, wind, 50.0)

# The wind farm state matrix is now updated with new TI values
current_ti = wf.States_WF[wf.StartI, 3]
```

# Notes
- The function modifies the wind farm object in-place (indicated by the `!` suffix)
- This "no correction" approach provides baseline turbulence intensity without 
  applying wake-induced corrections or measurement-based adjustments
- Error handling ensures graceful failure if matrix dimensions are incompatible

# See Also
- [`getDataTI`](@ref): Function used to retrieve turbulence intensity data
- [`TI_None`](@ref): Dispatch type for no correction strategy
"""
function correctTI!(::TI_None, set::Settings, wf::WindFarm, wind::Wind, t)
    # correctTI! updates the turbulent intensity (TI) value inT.States_WF
    # at the rowT.StartI and column 3, using data from getDataTI.

    # Get new turbulent intensity value
    TI = getDataTI(set, wind, wf, t)

    # Update the TI value in the States_WF matrix
    try
       wf.States_WF[wf.StartI, 3] .= TI'
    catch e
        error("Error updating wf.States_WF: $(e.msg)")
    end

    return nothing
end

"""
    getDataTI(set::Settings, wind::Wind, wf::WindFarm, t) -> Vector

Retrieve turbulence intensity data for all turbines at the current simulation time.

This function obtains turbulence intensity values for all turbines in the wind farm
using the configured turbulence model and wind field data. It serves as a wrapper
around the underlying turbulence intensity retrieval system.

# Arguments
- `set::Settings`: Settings object containing simulation configuration
  - `set.turb_mode`: Turbulence model configuration specifying the retrieval method
- `wind::Wind`: Wind configuration object containing turbulence intensity data
  - `wind.ti`: Turbulence intensity data, parameters, or model configuration
- `wf::WindFarm`: Wind farm object containing turbine information
  - `wf.nT`: Number of turbines in the wind farm
- `t`: Current simulation time for time-dependent turbulence intensity models

# Returns
- `Vector`: Turbulence intensity values for all turbines (dimensionless, typically 0.05-0.25)

# Behavior
The function creates a vector of all turbine indices `[1, 2, ..., nT]` and retrieves
the corresponding turbulence intensity values using the specified turbulence model.
The actual retrieval method depends on the `set.turb_mode` configuration and can include:
- Constant turbulence intensity
- Time-interpolated values from data files
- Turbine-specific interpolation
- Random walk models with covariance

# Example
```julia
# Get turbulence intensity for all turbines at t=100s
TI_values = getDataTI(settings, wind_config, wind_farm, 100.0)
println("TI for turbine 1: ", TI_values[1])
```

# See Also
- [`getWindTiT`](@ref): Underlying function for turbulence intensity retrieval
- [`correctTI!`](@ref): Function that uses this data to update wind farm states
"""
function getDataTI(set::Settings, wind::Wind, wf::WindFarm, t)
    TI = getWindTiT(set.turb_mode, wind.ti, collect(1:wf.nT), t)

    return TI
end

