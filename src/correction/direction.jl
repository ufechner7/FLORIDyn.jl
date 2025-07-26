# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function getDataDir(set::Settings, Wind, wf, SimTime)
    # Reads wind data and returns the current phi for all turbines

    if Wind.input_dir == "RW_with_Mean"
        phi = getWindDirT(set.dir_mode,wf.States_WF[wf.StartI, 2], Wind.dir)
    else
        phi = getWindDirT(set.dir_mode, Wind.dir, collect(1:wf.nT), SimTime)
    end

    return phi
end

"""
    correctDir!(::Direction_All, set::Settings, wf, Wind, t)

Corrects the direction based on the provided parameters.

# Arguments
- `::Direction_All`: The direction correction strategy or type.
- `set`:     The settings for the simulation.`
- `wf`:      The [WindFarm](@ref)
- `wind`:    The wind data or wind state.
- `t`:       The simulation time.

# Description
This function applies a direction correction using the specified strategy, updating the state in-place.
"""
function correctDir!(::Direction_All, set::Settings, wf, wind, t)
    # Get Data
    phi = getDataDir(set, wind, wf, t)
    # Correct
    wf.States_WF[:, 2] .= phi[1]
    # OP Orientation = turbine wind direction
    if size(wf.States_WF, 2) == 4
       wf.States_WF[wf.StartI, 4] .= phi[1]
    end
    return nothing
end

# function correctDir!(::Direction_All, set::Settings, wf, Wind, SimTime)
#     # CORRECTDIR Correction of the wind direction

#     ## Get Data
#     phi = getDataDir(set, Wind, wf, SimTime)

#     ## Correct the wind direction in the turbine states
#     wf.States_WF[:, 2] .= phi[1]

#     # OP Orientation = turbine wind direction
#     if size(wf.States_WF, 2) == 4
#        wf.States_WF[wf.StartI, 4] .= phi[1]
#     end

#     return wf
# end

