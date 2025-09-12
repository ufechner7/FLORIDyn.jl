# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Calculate axial induction factor, and calculate the demand

using FLORIDyn, ControlPlots, YAML

settings_file = get_default_project()[2]

# get the settings for the wind field, simulator, controller and turbine array
wind, sim, con, floris, floridyn, ta = setup(settings_file)


time_step = 4.0 # seconds

function calc_axial_induction(turbine_group, time)
    # Example: simple constant axial induction factor
    return 1/3  # Constant value for all turbines
end

function calc_axial_induction(turbine, time)
    turbine_group = turbine_group(ta, turbine)
    return 1/3  # Constant value for all turbines
end

function calc_demand(time)
    # Example: simple constant demand
    return 1.0  # Constant value
end