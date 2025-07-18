# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using FLORIDyn, YAML

settings_file = "data/2021_9T_Data.yaml"

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn = setup(settings_file)

# create settings struct
set = Settings(wind, sim)

# % Load linked data
turbProp        = turbineArrayProperties(settings_file)
paramFLORIS     = floris
paramFLORIDyn   = floridyn

T, wind, sim, con, paramFLORIS = prepareSimulation(wind, con, paramFLORIDyn, paramFLORIS, turbProp, sim)
