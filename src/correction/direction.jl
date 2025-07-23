# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function getDataDir(set::Settings, Wind, wf, SimTime)
    # Reads wind data and returns the current phi for all turbines

    if Wind.input_dir == "RW_with_Mean"
        phi = getWindDirT(set.dir_mode,wf.States_WF[wf.StartI, 2], Wind.dir)
    elseif Wind.input_dir == "EnKF_ZOH"
        phi =wf.States_WF[wf.StartI, 2]
    elseif Wind.input_dir == "EnKF_RW"
        phi =wf.States_WF[wf.StartI, 2]
        # (randn(1, length(phi))*Wind.dir.CholSig)' in MATLAB
        # randn generates standard normals. In Julia, randn(length) gives a vector.
        # Matrix multiplication and transpose need to be handled explicitly.
        phi = phi + (randn(length(phi))' * Wind.dir.CholSig)'
    elseif Wind.input_dir == "EnKF_InterpTurbine"
        # (1:wf.nT)' in MATLAB is just 1:wf.nT in Julia if used as a range, need to convert to an array if needed.
        phi = getWindDirT_EnKF(set.dir_mode, Wind.dir, collect(1:wf.nT), SimTime)
    elseif Wind.input_dir == "CLC_weighted_ZOH"
        # wf.C_Dir *wf.States_WF(:,2) in MATLAB is wf.C_Dir *wf.States_WF[:,2] in Julia
        phi = wf.C_Dir *wf.States_WF[:, 2]
    else
        phi = getWindDirT(set.dir_mode, Wind.dir, collect(1:wf.nT), SimTime)
    end

    return phi
end

"""
    correctDir!(::Direction_All, set::Settings, wf, Wind, SimTime)

Corrects the direction based on the provided parameters.

# Arguments
- `::Direction_All`: The direction correction strategy or type.
- `set`: The settings for the simulation.`
- `wf`: The current turbine (???)
- `Wind`: The wind data or wind state.
- `SimTime`: The simulation time.

# Description
This function applies a direction correction using the specified strategy, updating the state in-place.
"""
function correctDir!(::Direction_Interpolation, set::Settings, wf, Wind, SimTime)
    # Get Data
    phi = getDataDir(set, Wind, wf, SimTime)
    # Correct
   wf.States_WF[:, 2] .= phi[1]
    # OP Orientation = turbine wind direction
    if size(wf.States_WF, 2) == 4
       wf.States_WF[wf.StartI, 4] .= phi[1]
    end
    return nothing
end

function correctDir!(::Direction_All, set::Settings, wf, Wind, SimTime)
    # CORRECTDIR Correction of the wind direction

    ## Get Data
    phi = getDataDir(set, Wind, wf, SimTime)

    ## Correct the wind direction in the turbine states
   wf.States_WF[:, 2] .= phi[1]

    # OP Orientation = turbine wind direction
    if size(wf.States_WF, 2) == 4
       wf.States_WF[wf.StartI, 4] .= phi[1]
    end

    return wf
end

