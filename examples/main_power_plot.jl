# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# MainFLORIDyn Center-Line model
# Improved FLORIDyn approach over the gaussian FLORIDyn model

# Minimal example of how to run a simulation using FLORIDyn.jl
using FLORIDyn, TerminalPager, DistributedNext, Statistics
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"
PLOT_FLOW_FIELD = false

# Load vis settings from YAML file
vis = Vis(vis_file)
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end
pltctrl = nothing
# Provide ControlPlots module only for pure sequential plotting (single-threaded, no workers)
if Threads.nthreads() == 1
    pltctrl = ControlPlots
end

# Automatic parallel/threading setup
include("../examples/remote_plotting.jl")

times, rel_power, set, wf, wind, floris = calc_rel_power(settings_file; dt=350, wind_dir=nothing)

plot_rmt(times, rel_power .* 100; xlabel="Time [s]", ylabel="Rel. Power Output [%]", pltctrl)

println("\nMean Relative Power Output:  $(round((mean(rel_power) * 100), digits=2)) %")
println("Final Relative Power Output: $(round((rel_power[end] * 100), digits=2)) %")

if PLOT_FLOW_FIELD
    Z, X, Y = calcFlowField(set, wf, wind, floris; plt, vis)
    plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)
end
