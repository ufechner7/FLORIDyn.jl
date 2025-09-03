# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# This script shows how to plot multiple time series in one plot in the remote process (rmt)

using LaTeXStrings, DistributedNext, FLORIDyn
if Threads.nthreads() == 1; using ControlPlots; end

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

X  = 0:0.1:2pi
Y1 = sin.(X)
Y2 = cos.(X)
plot_rmt(X, [Y1, Y2]; xlabel=L"\alpha = [0..2\pi]", ylabel="Function Value", labels=["sin","cos"], xlims=(2, 6), 
         fig="Dual", pltctrl)