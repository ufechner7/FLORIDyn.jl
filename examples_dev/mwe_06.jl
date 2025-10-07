# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example how to plot two time series in one plot with one y axis using plot_rmt().
using FLORIDyn, DistributedNext
using ControlPlots

# Select plotting backend: use ControlPlots only in single-threaded mode; otherwise let plot_rmt handle remote plotting
pltctrl = nothing
if Threads.nthreads() == 1
    pltctrl = ControlPlots
end

# Initialize remote plotting when multithreaded so plot_rmt can dispatch to a worker process
include("../examples/remote_plotting.jl")

# Create example data
time = 1:10
series1 = rand(10)
series2 = rand(10)

plot_rmt(time, [series1, series2]; 
         xlabel="Time", 
         ylabel="Values",
         labels=["Series 1", "Series 2"], 
         title="Two Time Series",
         pltctrl)