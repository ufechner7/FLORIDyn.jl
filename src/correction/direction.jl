"""
    correctDir!(::Direction_All, T, Wind, SimTime)

Corrects the direction based on the provided parameters.

# Arguments
- `::Direction_All`: The direction correction strategy or type.
- `T`: The current turbine (???)
- `Wind`: The wind data or wind state.
- `SimTime`: The simulation time.

# Description
This function applies a direction correction using the specified strategy, updating the state in-place.
"""
function correctDir!(::Direction_All, T, Wind, SimTime)
    # Get Data
    phi = getDataDir(Wind, T, SimTime)
    # Correct
    T.States_WF[:, 2] .= phi[1]
    # OP Orientation = turbine wind direction
    if size(T.States_WF, 2) == 4
        T.States_WF[T.StartI, 4] = phi[1]
    end
    return nothing
end
