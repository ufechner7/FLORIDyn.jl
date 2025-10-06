# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Calculate axial induction factor, and calculate the demand

using FLORIDyn, ControlPlots, YAML, DistributedNext

# Select plotting backend: use ControlPlots only in single-threaded mode; otherwise let plot_rmt handle remote plotting
pltctrl = nothing
if Threads.nthreads() == 1
    pltctrl = ControlPlots
end

# Initialize remote plotting when multithreaded so plot_rmt can dispatch to a worker process
include("remote_plotting.jl")
include("calc_induction_matrix.jl")

USE_TGC = true
USE_FEED_FORWARD = true
USE_STEP = true

settings_file = get_default_project()[2]

# get the settings for the wind field, simulator, controller and turbine array
wind, sim, con, floris, floridyn, ta = setup(settings_file)
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Calculate time step and relative end time based on simulation settings
time_step = sim.time_step  # seconds
t_end = sim.end_time - sim.start_time  # relative end time in seconds

function plot_demand()
    # Create time vector from 0 to t_end with time_step intervals
    time_vector = 0:time_step:t_end

    # Calculate demand for each time point
    demand_values = [calc_demand(t) for t in time_vector]
    
    # Calculate actual induction values for representative turbines from each group
    # (This includes the time-dependent corrections implemented in calc_axial_induction)
    
    # Find representative turbines for each group
    turbine_group1 = 1   # Turbine 1 is in group 1
    turbine_group2 = 2   # Turbine 2 is in group 2
    turbine_group3 = 3   # Turbine 3 is in group 3
    turbine_group4 = 7   # Turbine 7 is in group 4
    
    # Calculate induction values using actual calc_axial_induction (with corrections)
    induction_group1 = [calc_axial_induction(ta, con, turbine_group1, t) for t in time_vector]
    induction_group2 = [calc_axial_induction(ta, con, turbine_group2, t) for t in time_vector]
    induction_group3 = [calc_axial_induction(ta, con, turbine_group3, t) for t in time_vector]
    induction_group4 = [calc_axial_induction(ta, con, turbine_group4, t) for t in time_vector]
    
    # Combine all data series
    all_data = [demand_values, induction_group1, induction_group2, induction_group3, induction_group4]
    labels = ["Demand", "Group 1", "Group 2", "Group 3", "Group 4"]
    
    # Create the plot (thread-safe via plot_rmt)
    plot_rmt(time_vector, all_data;
             xlabel="Time [s]",
             ylabel="Demand and Induction [-]",
             title="Demand and Induction vs Time for All Groups",
             labels=labels,
             pltctrl=pltctrl)
end

function plot_induction_matrix()
    # Calculate induction matrix for all turbines over time
    induction_matrix = calc_induction_matrix(ta, con, time_step, t_end)
    time_vector = induction_matrix[:, 1]  # Extract time from first column
    n_time_steps = size(induction_matrix, 1)
    
    # Initialize group data arrays
    group_data = [Float64[] for _ in 1:4]  # Arrays for groups 1-4
    
    # Calculate average induction for each group at each time step
    for t_idx in 1:n_time_steps
        group_sums = zeros(4)
        group_counts = zeros(Int, 4)
        
        # Sum induction values for each group (columns 2 onwards are turbine data)
        for turbine in 1:size(ta.pos, 1)
            group_id = FLORIDyn.turbine_group(ta, turbine)
            if 1 <= group_id <= 4
                group_sums[group_id] += induction_matrix[t_idx, turbine + 1]  # +1 because first column is time
                group_counts[group_id] += 1
            end
        end
        
        # Calculate averages and store
        for group in 1:4
            if group_counts[group] > 0
                avg_induction = group_sums[group] / group_counts[group]
            else
                avg_induction = 0.0
            end
            push!(group_data[group], avg_induction)
        end
    end
    
    # Create labels with group information
    group_labels = ["Group 1", "Group 2", "Group 3", "Group 4"]
    
    plot_rmt(time_vector, group_data;
             xlabel="Time [s]",
             ylabel="Axial Induction Factor [-]",
             title="Average Axial Induction Factor vs Time by Turbine Group",
             labels=group_labels,
             pltctrl=pltctrl)
end

con.induction_data = calc_induction_matrix(ta, con, time_step, t_end)

plot_induction_matrix()
