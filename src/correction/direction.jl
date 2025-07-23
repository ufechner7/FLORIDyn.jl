# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function getDataDir(set::Settings, Wind, T, SimTime)
    # Reads wind data and returns the current phi for all turbines

    if Wind.input_dir == "RW_with_Mean"
        phi = getWindDirT(set.dir_mode,T.States_WF[T.StartI, 2], Wind.dir)
    elseif Wind.input_dir == "EnKF_ZOH"
        phi =T.States_WF[T.StartI, 2]
    elseif Wind.input_dir == "EnKF_RW"
        phi =T.States_WF[T.StartI, 2]
        # (randn(1, length(phi))*Wind.dir.CholSig)' in MATLAB
        # randn generates standard normals. In Julia, randn(length) gives a vector.
        # Matrix multiplication and transpose need to be handled explicitly.
        phi = phi + (randn(length(phi))' * Wind.dir.CholSig)'
    elseif Wind.input_dir == "EnKF_InterpTurbine"
        # (1:T.nT)' in MATLAB is just 1:T.nT in Julia if used as a range, need to convert to an array if needed.
        phi = getWindDirT_EnKF(set.dir_mode, Wind.dir, collect(1:T.nT), SimTime)
    elseif Wind.input_dir == "CLC_weighted_ZOH"
        # T.C_Dir *T.States_WF(:,2) in MATLAB is T.C_Dir *T.States_WF[:,2] in Julia
        phi = T.C_Dir *T.States_WF[:, 2]
    else
        phi = getWindDirT(set.dir_mode, Wind.dir, collect(1:T.nT), SimTime)
    end

    return phi
end

"""
    correctDir!(::Direction_All, set::Settings, T, Wind, SimTime)

Corrects the direction based on the provided parameters.

# Arguments
- `::Direction_All`: The direction correction strategy or type.
- `set`: The settings for the simulation.`
- `T`: The current turbine (???)
- `Wind`: The wind data or wind state.
- `SimTime`: The simulation time.

# Description
This function applies a direction correction using the specified strategy, updating the state in-place.
"""
function correctDir!(::Direction_Interpolation, set::Settings, T, Wind, SimTime)
    # Get Data
    phi = getDataDir(set, Wind, T, SimTime)
    # Correct
   T.States_WF[:, 2] .= phi[1]
    # OP Orientation = turbine wind direction
    if size(T.States_WF, 2) == 4
       T.States_WF[T.StartI, 4] .= phi[1]
    end
    return nothing
end

function correctDir!(::Direction_All, set::Settings, T, Wind, SimTime)
    # CORRECTDIR Correction of the wind direction

    ## Get Data
    phi = getDataDir(set, Wind, T, SimTime)

    ## Correct the wind direction in the turbine states
   T.States_WF[:, 2] .= phi[1]

    # OP Orientation = turbine wind direction
    if size(T.States_WF, 2) == 4
       T.States_WF[T.StartI, 4] .= phi[1]
    end

    return T
end

