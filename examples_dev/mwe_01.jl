# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Check if the function getMeasurements is type stable.
# The output must NOT contain any RED variables.
# For performance reasons, the function should be type stable.

using FLORIDyn

function get_parameters()
    settings_file = "data/2021_9T_Data.yaml"

    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)

    # create settings struct
    set = Settings(wind, sim, con)

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)  
    wf = initSimulation(wf, sim)
    
    # Create visualization settings for testing
    vis = Vis(online=false, save=false, unit_test=true)
    wf, md, mi = runFLORIDyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)
    return wf, set, floris, wind, md
end

# Create minimal test matrices
mx = [0.0 100.0; 200.0 300.0]  # 2x2 grid of x-coordinates
my = [0.0 0.0; 100.0 100.0]    # 2x2 grid of y-coordinates
nM = 3  # Number of measurements (typically velocity reduction, added turbulence, effective wind speed)
zh = 90.0  # Hub height
wf, set, floris, wind, md = get_parameters()
# Call the function with the correct number of parameters
@code_warntype getMeasurements(mx, my, nM, zh, wf, set, floris, wind)

nothing