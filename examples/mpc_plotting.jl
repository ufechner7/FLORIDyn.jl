# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# This file is included by mpc.jl and contains the plotting functions.

"""
    plot_induction(vis, optimal_correction::Vector{Float64})

Plot the axial induction factor for turbine group 1 over time range 500-1500s.

# Arguments
- `vis`: Visualization object containing `t_skip` parameter
- `optimal_correction::Vector{Float64}`: Optimal correction parameters from optimization

# Description
Uses `calc_axial_induction2` to compute induction values for group 1 and plots
them against time. The plot uses the global `pltctrl` variable for thread-safe plotting.
"""
function plot_induction(vis, optimal_correction::Vector{Float64})
    # Time range: 500 to 1500 seconds
    t_start = vis.t_skip
    t_end   = vis.t_skip + T_END + T_EXTRA
    dt = time_step
    
    # Print diagnostic information
    println("\n=== Diagnostic: plot_induction ===")
    println("optimal_correction[1:$CONTROL_POINTS]: ", optimal_correction[1:CONTROL_POINTS])
    if length(optimal_correction) > CONTROL_POINTS
        println("optimal_correction[$(CONTROL_POINTS+1)] (group 1 id_correction): ", optimal_correction[CONTROL_POINTS+1])
    end
    
    # Create time vector
    time_vec = t_start:dt:t_end
    n_points = length(time_vec)
    
    # Calculate induction for group 1 at each time point
    induction_values = zeros(n_points)
    group_id = 1
    
    for (i, t) in enumerate(time_vec)
        induction, _ = calc_axial_induction2(vis, t, optimal_correction; group_id=group_id)
        induction_values[i] = induction
    end
    
    # Find and report the dip
    min_idx = argmin(induction_values)
    min_time = collect(time_vec)[min_idx]
    min_val = induction_values[min_idx]
    println("Minimum induction: $(round(min_val, digits=5)) at time $(round(min_time, digits=1))s")
    println("==================================\n")
    
    # Plot
    plot_rmt(collect(time_vec), induction_values;
             xlabel="Time [s]",
             ylabel="Axial Induction Factor [-]",
             title="Induction Factor for Turbine Group 1",
             fig="Induction Group 1",
             pltctrl=pltctrl)
end

"""
    plot_correction_curve(optimal_correction::Vector{Float64})

Plot the correction curve from the piecewise cubic Hermite spline interpolation over s=0..1.

# Arguments
- `optimal_correction::Vector{Float64}`: Optimal correction parameters from optimization

# Description
Plots the piecewise cubic Hermite spline interpolation curve showing how the correction
factor varies across the normalized parameter s from 0 to 1. Uses the first CONTROL_POINTS
elements of `optimal_correction` as control points evenly spaced from s = 0 to s = 1.0.
"""
function plot_correction_curve(optimal_correction::Vector{Float64})
    # Create s vector from 0 to 1
    s_vec = 0.0:0.01:1.0
    n_points = length(s_vec)
    
    # Calculate correction_result for each s value
    correction_values = zeros(n_points)
    
    for (i, s) in enumerate(s_vec)
        correction_values[i] = interpolate_hermite_spline(s, optimal_correction[1:CONTROL_POINTS])
    end
    
    # Print diagnostic information
    println("\n=== Diagnostic: plot_correction_curve ===")
    println("Control points (correction[1:$CONTROL_POINTS]): ", optimal_correction[1:CONTROL_POINTS])
    println("Min correction: $(round(minimum(correction_values), digits=4))")
    println("Max correction: $(round(maximum(correction_values), digits=4))")
    println("======================================\n")
    
    # Plot
    plot_rmt(collect(s_vec), correction_values;
             xlabel="Normalized Parameter s [-]",
             ylabel="Correction Factor [-]",
             title="Hermite Spline Interpolation Correction Curve",
             fig="Correction Curve",
             pltctrl=pltctrl)
end
function plot_correction2(optimal_correction::Vector{Float64})
    plt.figure("Hermite Spline Interpolation Correction and Control Points")
    points = optimal_correction[1:CONTROL_POINTS]
    s_vec=0:CONTROL_POINTS-1
    plt.plot(s_vec./(CONTROL_POINTS-1), points, label="Control Points",marker="x", linestyle="None")
    # Create s vector from 0 to 1
    s_vec1 = 0.0:0.01:1.0
    n_points = length(s_vec1)
    
    # Calculate correction_result for each s value
    correction_values = zeros(n_points)
    
    for (i, s) in enumerate(s_vec1)
        correction_values[i] = interpolate_hermite_spline(s, optimal_correction[1:CONTROL_POINTS])
    end
    plt.plot(collect(s_vec1), correction_values, label="Interpolated", linestyle="--")
    plt.xlabel("Normalized Parameter s [-]")
    plt.ylabel("Correction Factor [-]")
    plt.title("Hermite Spline Interpolation Correction and Control Points")
    plt.legend()
    plt.grid(true, color = "#DDDDDD")
    plt.show()
end

function plot_power_and_demand(time_vector, rel_power, rel_power_ref, demand_values; vis, pltctrl)
    if GROUP_CONTROL
        plot_rmt(time_vector, [rel_power[1:length(time_vector)] .* 100, rel_power_ref[1:length(time_vector)] .* 100, demand_values .* 100]; xlabel="Time [s]", xlims=(vis.t_skip, time_vector[end]),
                ylabel="Rel. Power Output [%]", labels=["rel_power", "rel_power_ref", "rel_demand"], title="Rel. Power and Demand "*string(GROUPS)*" TGs", fig="Rel. Power and Demand", pltctrl)
    else
        plot_rmt(time_vector, [rel_power[1:length(time_vector)] .* 100, demand_values .* 100]; xlabel="Time [s]", xlims=(vis.t_skip, time_vector[end]),
                ylabel="Rel. Power Output [%]", labels=["rel_power", "rel_demand"], fig="Rel. Power and Demand", title="Rel. Power and Demand (1 TG)", pltctrl)
    end
end

function plot_axial_induction()
    # Plot average axial induction factor per turbine group over time
    # Extract time from first column of induction_data
    time_vec_ind = induction_data[:, 1]
    n_time_steps = size(induction_data, 1)
    n_turbines = size(ta.pos, 1)

    # Prepare containers for the groups
    group_data = [Float64[] for _ in 1:GROUPS]

    # Since all turbines in a group have identical induction, pick one representative turbine per group
    group_indices = [findfirst(i -> FLORIDyn.turbine_group(ta, i) == g, 1:n_turbines) for g in 1:GROUPS]
    for g in 1:GROUPS
        idx = group_indices[g]
        if isnothing(idx)
            group_data[g] = fill(0.0, n_time_steps)
        else
            group_data[g] = induction_data[:, idx + 1]
        end
    end
    
    # Create group labels dynamically
    group_labels = ["Group $i" for i in 1:GROUPS]
    
    plot_rmt(time_vec_ind, group_data;
             xlabel="Time [s]",
             ylabel="Axial Induction Factor [-]",
             title="Average Axial Induction Factor vs Time by Turbine Group",
             labels=group_labels,
             fig="Induction by Group",
             pltctrl=pltctrl)
    
    # Create bar plot of average induction factor per group
    avg_induction = [isempty(group_data[g]) ? 0.0 : mean(group_data[g]) for g in 1:GROUPS]
    
    if !isnothing(plt)
        plt.figure(figsize=(10, 6))
        
        # Create color palette - cycle through colors if more than 8 groups
        base_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]
        colors = [base_colors[mod1(i, length(base_colors))] for i in 1:GROUPS]
        
        bars = plt.bar(1:GROUPS, avg_induction, color=colors)
        plt.xlabel("Turbine Group", fontsize=12)
        plt.ylabel("Average Axial Induction Factor [-]", fontsize=12)
        plt.title("Average Axial Induction Factor by Turbine Group", fontsize=14)
        plt.xticks(1:GROUPS, ["Group $i" for i in 1:GROUPS], rotation=45, ha="right")
        plt.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        
        # Add value labels on top of bars
        for (i, v) in enumerate(avg_induction)
            plt.text(i, v, @sprintf("%.3f", v), ha="center", va="bottom", fontsize=10)
        end
    end
end