# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model
using FLORIDyn, TerminalPager

settings_file = "data/2021_9T_Data.yaml"

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn = setup(settings_file)

# create settings struct
set = Settings(wind, sim, con)

# % Load linked data
turbProp        = turbineArrayProperties(settings_file)

T, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbProp, sim)

# Run initial conditions until no more change happens (wrong comment in original code)
T = initSimulation(T, wind, sim, con, floridyn, floris)

@time T, M, Mint = FLORIDynCL(set, T, wind, sim, con, floridyn, floris)
# 0.16 s on Desktop, 0.40 s with MATLAB

@info "Type 'M |> pager' to see the results of the simulation."
@info "Type 'q' to exit the pager."
nothing
