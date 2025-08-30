# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Minimal example of how to run a simulation using FLORIDyn.jl
using FLORIDyn, TerminalPager, DistributedNext, BenchmarkTools

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"

# Load vis settings from YAML file
vis = Vis(vis_file)
vis.show_plots = false
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Run initial conditions
wf = initSimulation(wf, sim)

vis.online = false
"""
Benchmark for a large 54 turbine case.

Measures ONLY the simulation loop (excludes setup, initialisation, compilation).
"""

# Warm-up (compile) — not timed by BenchmarkTools below
run_floridyn(plt, set, deepcopy(wf), wind, deepcopy(sim), con, vis, floridyn, floris; collect_md=false, collect_interactions=false)
GC.gc()

const initial_wf = deepcopy(wf)
const initial_sim = deepcopy(sim)

@btime begin
    run_floridyn($plt, $set, wf_copy, $wind, sim_copy, $con, $vis, $floridyn, $floris; collect_md=false, collect_interactions=false)
end setup=(wf_copy=deepcopy(initial_wf); sim_copy=deepcopy(initial_sim)) evals=1

nothing

# Previous: ~1.064 s including measurement allocations
# New benchmark isolates pure simulation runtime.
