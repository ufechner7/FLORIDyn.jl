# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example of a simple Julia script to demonstrate local plotting.

using Timers, ControlPlots, FLORIDyn

tic()
include("../examples/remote_plotting.jl") 
init_plotting()
toc()

X = 0:0.1:10
Y = sin.(X)

@time display(plot(X, Y, fig="Sine Function Plot", xlabel="X-axis", ylabel="Y-axis"))

# Results are now local-only.

nothing
