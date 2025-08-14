# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Testcase for bug https://github.com/ufechner7/FLORIDyn.jl/issues/35
using FLORIDyn, TerminalPager, MAT, ControlPlots, DistributedNext, LinearAlgebra, Statistics, DataFrames

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"
after_init_file = "test/data/after_init_simulation_T.mat"
T_init = matread(after_init_file)["T"]
turbines_init = (turbines(T_init))  # Creates a 1800×3 DataFrame with turbine states after init

function rel_err(a, b)
    return norm(a - b) / norm(b)
end

# Load vis settings from YAML file
vis = Vis(vis_file)

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
sim.n_sim_steps = 195

# compare wf.windfield with Matlab before running floridyn

turbines_wf = wf.turbines
function plot_dfs(df1, df2)
    p = plot(1:length(df1.yaw), [df1.TI, df2.TI]; 
             ylabel="TI", labels=["wf (Julia)", "Init"], fig="turbines(wf), turbines(Init)")
    display(p)
end
plot_dfs(turbines_init, turbines_wf)

# println("Relative error (turbines): ", round(rel_err(turbines_wf, turbines_ref)*100, digits=2), " %")
df1, df2 = compare_dataframes(turbines_wf, turbines_init)
println("Number of differing rows found: ", size(df1, 1), " out of ", size(turbines_wf, 1))


nothing
