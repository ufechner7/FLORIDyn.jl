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

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbProp, sim)

# Run initial conditions until no more change happens (wrong comment in original code)
wf = initSimulation(wf, sim)

@time wf, md, mi = runFLORIDyn(set, wf, wind, sim, con, floridyn, floris)
# 0.115 s on Desktop, 0.39 s with MATLAB
# 0.115 seconds (891.24 k allocations: 368.147 MiB, 10.18% gc time)
# 0.110 seconds (883.14 k allocations: 272.819 MiB, 9.62% gc time) iterateOPs! allocation free
# 0.081 seconds (723.31 k allocations: 168.226 MiB, 8.10% gc time) findTurbineGroups allocation free

@info "Type 'md |> pager' to see the results of the simulation."
@info "Type 'q' to exit the pager."
nothing
