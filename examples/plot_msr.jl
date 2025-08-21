# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, TOML, DistributedNext
if Threads.nthreads() == 1; 
    using ControlPlots
    v = VersionNumber(TOML.parsefile(joinpath(Base.pkgdir(ControlPlots), "Project.toml"))["version"])
    @assert v >= v"0.2.7" "This script requires ControlPlots version 0.2.7 or higher."
end

MULTI = true

settings_file, vis_file = get_default_project()

vis = Vis(vis_file)
if (@isdefined plt) && !isnothing(plt)
    plt.ion()
else
    plt = nothing
end
plt1 = nothing
if Threads.nthreads() == 1; plt1 = ControlPlots; end

# Automatic parallel/threading setup
include("remote_plotting.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# create settings struct with automatic threading/parallel detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)

wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
vis.online = false
@time wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

# Arrays to store time series data
times = Float64[]
measurements = Vector{Float64}[]
turbines = 1:wf.nT

msr =  get_default_msr()

# Determine measurement type based on msr parameter
if msr == VelReduction
    data_column = "ForeignReduction"
    title = "Velocity Reduction"
    ylabel = "Vel Reduction [%]"
    msr_name = "msr_velocity_reduction"
elseif msr == AddedTurbulence
    data_column = "AddedTurbulence" 
    title = "Added Turbulence"
    ylabel = "Added Turbulence [%]"
    msr_name = "msr_added_turbulence"
elseif msr == EffWind
    data_column = "EffWindSpeed"
    title = "Effective Wind Speed"
    ylabel = "Effective Wind Speed [m/s]"
    msr_name = "msr_eff_wind_speed"
else
    error("Invalid msr value: $msr. Must be VelReduction, AddedTurbulence, or EffWind.")
end

# Check if the required column exists in the DataFrame
if !(data_column in names(md))
    error("Column '$data_column' not found in measurement data. Available columns: $(names(md))")
end
measurement_data = md[!, data_column]
timeFDyn = md.Time .- md.Time[1]
lay = get_layout(wf.nT)

for iT in 1:wf.nT
    push!(measurements,  measurement_data[iT:wf.nT:end])
end

for time in sim.start_time:sim.time_step:sim.end_time
    push!(times, time - sim.start_time)
end

# # Convert vector of vectors to matrix for easier plotting
msr_matrix = hcat(measurements...)  # Transpose to get time × turbines
msr_matrix = msr_matrix'          # Now it's time × turbines

# # Create dynamic plot arguments based on number of turbines
n_turbines = wf.nT


rows, lines = get_layout(wf.nT)

# Group turbines into subplots based on layout
plot_data = []
turbine_labels = []
subplot_labels = []
    
turbine_idx = 1
for row in 1:rows
    global turbine_idx, lines_in_subplot, labels_in_subplot
    if turbine_idx > n_turbines
        break
    end
        
    # Collect lines for this subplot
    lines_in_subplot = Vector{Vector{Float64}}()
    labels_in_subplot = Vector{String}()

    for line in 1:lines
        if turbine_idx <= n_turbines
            push!(lines_in_subplot, msr_matrix[:, turbine_idx])
            push!(labels_in_subplot, "T$(turbines[turbine_idx])")
            turbine_idx += 1
        end
    end
        
    # Add subplot data
    if length(lines_in_subplot) == 1
        push!(plot_data, lines_in_subplot[1])
    else
        push!(plot_data, lines_in_subplot)
    end
    
    push!(turbine_labels, "Wind Direction [°]")
    push!(subplot_labels, labels_in_subplot)
end
    
# Plot with multiple lines per subplot
plot_x(times, plot_data...; ylabels=turbine_labels, labels=subplot_labels,
            fig=title, xlabel="rel_time [s]", ysize = 10, bottom=0.02, plt=plt1)


nothing