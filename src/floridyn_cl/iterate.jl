# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    iterateOPs!(::IterateOPs_basic, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn)

Advance operational points through the wind field using time-marching dynamics.

This function implements the basic time-stepping algorithm for operational points (OPs) in the 
FLORIDyn simulation. It advances all operational points through space based on wind velocities 
and wake deflection effects, maintaining proper ordering and updating turbine states accordingly.

# Arguments
- `::IterateOPs_basic`: Dispatch type indicating the basic iteration algorithm
- `wf::WindFarm`: Wind farm object containing all turbine and operational point states. See [`WindFarm`](@ref)
  - `wf.States_OP`: Matrix of operational point states [x, y, z, `dw_pos`, `cw_x`, `cw_z`, ...]
  - `wf.States_T`: Matrix of turbine states 
  - `wf.States_WF`: Matrix of wind field states [velocity, direction, `turbulence_intensity`, ...]
  - `wf.StartI`: Starting indices for each turbine's operational points
  - `wf.nT`: Number of turbines
  - `wf.nOP`: Number of operational points per turbine
  - `wf.D`: Rotor diameters [m]
- `sim::Sim`: Simulation configuration object. See [`Sim`](@ref)
  - `sim.time_step`: Time step size [s]
  - `sim.dyn.advection`: Advection scaling factor
- `floris::Floris`: FLORIS model parameters for wake deflection calculations. See [`Floris`](@ref)
- `floridyn::FloriDyn`: FLORIDyn model parameters. See [`FloriDyn`](@ref)

# Returns
- nothing

# Algorithm
The function performs the following steps for each time iteration:

## 1. State Preservation
Saves the current turbine operational point states before advancement to preserve boundary conditions.

## 2. Spatial Advancement
### Downwind Movement
Advances operational points downstream based on local wind velocity:
```
step_dw = Δt × U × advection_factor
OP_dw_position += step_dw
```

### Crosswind Deflection
Calculates wake centerline deflection using FLORIS model and updates crosswind positions:
```
deflection = centerline(States_OP, States_T, States_WF, floris, D)
OP_cw_position = deflection
```

## 3. Coordinate System Transformation
Converts local wind-aligned movements to world coordinates using wind direction:
```
φ = angSOWFA2world(wind_direction)
x_world += cos(φ) × step_dw - sin(φ) × step_cw_x
y_world += sin(φ) × step_dw + cos(φ) × step_cw_x
z_world += step_cw_z
```

## 4. Temporal Shifting
Uses circular shifting to advance the time history:
- Shifts all state matrices by one time step
- Initializes new operational points with saved turbine states
- Maintains temporal continuity of the simulation

## 5. Spatial Ordering
Ensures operational points remain ordered by downstream position for each turbine:
- Sorts operational points by downstream position (`States_OP[:, 4]`)
- Maintains consistency across all state matrices

# Mathematical Description
The coordinate transformation from wind-aligned to world coordinates follows:
```
[x']   [cos(φ)  -sin(φ)] [step_dw ]
[y'] = [sin(φ)   cos(φ)] [step_cw_x]
```
where `φ` is the wind direction angle in world coordinates.

# Notes
- The function modifies the wind farm object in-place (indicated by the `!` suffix)
- Temporal shifting maintains a moving window of operational point history
- Spatial ordering ensures downstream distance monotonicity for wake calculations
- The algorithm handles both advection and deflection physics simultaneously
- Coordinate transformations account for SOWFA wind direction conventions via [`angSOWFA2world`](@ref)
"""
@views function iterateOPs!(::IterateOPs_basic, wf::WindFarm, sim::Sim, floris::Floris, floridyn::FloriDyn)
    # Save turbine OPs
    tmpOPStates = copy(wf.States_OP[wf.StartI, :])
    tmpTStates  = copy(wf.States_T[wf.StartI, :])
    tmpWFSTates = copy(wf.States_WF[wf.StartI, :])

    # Shift states
    # Downwind step
    step_dw = sim.time_step .* sim.dyn.advection .* wf.States_WF[:, 1] 
     wf.States_OP[:, 4] .+= step_dw

    # Crosswind step
    deflection = centerline(wf.States_OP, wf.States_T, wf.States_WF, floris, wf.D[1])
    step_cw = deflection .- wf.States_OP[:, 5:6]
     wf.States_OP[:, 5:6] .= deflection

    # World coordinate system adjustment
    @inbounds @simd for i in axes(wf.States_OP,1)
      phiwi = angSOWFA2world(wf.States_WF[i, 2])  # Convert wind direction to world coordinates
      cphiwi = cos(phiwi)
      sphiwi = sin(phiwi)
      ai = step_cw[i, 1]
      sdwi = step_dw[i]
      wf.States_OP[i, 1] += cphiwi * sdwi - sphiwi * ai
      wf.States_OP[i, 2] += sphiwi * sdwi + cphiwi * ai
      wf.States_OP[i, 3] += step_cw[i, 2]
    end

    # Circshift & init first OPs
    # OPs
    wf.States_OP = circshift(wf.States_OP, (1, 0))
     wf.States_OP[wf.StartI, :] .= tmpOPStates

    # Turbines
    wf.States_T = circshift(wf.States_T, (1, 0))
     wf.States_T[wf.StartI, :] .= tmpTStates

    # Wind Farm
    wf.States_WF = circshift(wf.States_WF, (1, 0))
     wf.States_WF[wf.StartI, :] .= tmpWFSTates

    # Check if OPs are in order
    buf = zeros(Int,size(wf.States_OP, 1))
    for iT in 1:wf.nT
        inds = wf.StartI[iT]:(wf.StartI[iT] + wf.nOP - 1)
        indOP = buf[inds]
        sortperm!(indOP,wf.States_OP[inds, 4])
        if ! issorted(indOP)  # check if already sorted
           wf.States_OP[inds, :] .= wf.States_OP[inds[indOP], :]
           wf.States_T[inds, :]  .= wf.States_T[inds[indOP], :]
           wf.States_WF[inds, :] .= wf.States_WF[inds[indOP], :]
        end

    end
    return nothing
end
