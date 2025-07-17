# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using FLORIDyn, YAML

# get the settings for the wind field, simulator and controller
wind, sim, con = setup("data/2021_9T_Data.yaml")

# create settings struct
set = Settings(wind, sim)

# % Load linked data
# turbProp        = turbineArrayProperties();
