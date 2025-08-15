# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Testcase for bug https://github.com/ufechner7/FLORIDyn.jl/issues/35
using FLORIDyn, TerminalPager, MAT, ControlPlots, DistributedNext, LinearAlgebra, Statistics, DataFrames

settings_file = "data/2021_9T_Data.yaml"
vis_file      = "data/vis_default.yaml"
input_file    = "test/data/input_T_195_steps.mat"
T_ref = matread(input_file)["T"]
turbines_ref = (turbines(T_ref))  # Creates a 1800×3 DataFrame with turbine states
matlab_file   = "test/data/flowfield_xyz_195_steps.mat"
vars = matread(matlab_file)
X_ref = vars["X"]
Y_ref = vars["Y"]
Z_ref = vars["Z"]

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

vis.online = false
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

turbines_wf = wf.turbines

function plot_dfs(df1, df2; fig=nothing, max_op=195)
    # Filter dataframes to keep only rows where OP <= max_op
    df1_filtered = filter(row -> row.OP <= max_op, df1)
    df2_filtered = filter(row -> row.OP <= max_op, df2)
    
    turbine_dfs1 = [df for df in groupby(df1_filtered, :Turbine)]
    turbine_dfs2 = [df for df in groupby(df2_filtered, :Turbine)]
    
    # Create vectors for each turbine plot (each containing two lines: Julia and Ref)
    turbine_plots = [[collect(turbine_dfs1[i].TI), collect(turbine_dfs2[i].TI)] for i in 1:9]
    ylabels = ["TI Turbine $i" for i in 1:9]
    labels = [["Julia", "Ref"] for _ in 1:9]

    if isnothing(fig)
        fig = "Turbine TI Comparison"
    end
    
    p = plotx(collect(turbine_dfs1[1].OP), turbine_plots...; 
              xlabel="Operating Point", ylabels=ylabels, labels=labels, ysize=10,
              fig, bottom=0.02)
    display(p)
end

plot_dfs(turbines_wf, turbines_ref; fig="Turbine TI Comparison $(sim.n_sim_steps)")

max_op = 1  # Maximum operating point to consider   
turbines_wf = filter(row -> row.OP <= max_op, turbines_wf)
turbines_ref = filter(row -> row.OP <= max_op, turbines_ref)

df1, df2 = compare_dataframes(turbines_wf, turbines_ref)
println("Number of differing rows found: ", size(df1, 1), " out of ", size(turbines_wf, 1))


nothing
