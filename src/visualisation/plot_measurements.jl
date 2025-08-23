# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    plotMeasurements(plt, wf::WindFarm, md::DataFrame, vis::Vis; separated=false, msr=VelReduction) -> Nothing

Plot foreign reduction measurements from FLORIDyn simulation data.

# Arguments
- `plt`: Plotting package (e.g., PyPlot, which is exported from ControlPlots)
- `wf::WindFarm`: Wind farm object with field `nT` (number of turbines). See [`WindFarm`](@ref)
- `md::DataFrame`: Measurements DataFrame containing time series data with columns:
  - `Time`: Simulation time [s]
  - `ForeignReduction`: Foreign reduction  \\[%\\]
- `vis::Vis`: Visualization settings including unit_test parameter. See [`Vis`](@ref)
- `separated::Bool`: Whether to use separated subplot layout (default: false)
- `msr::MSR`: Measurement type to plot, see: [MSR](@ref) 

# Returns
- `nothing`

# Description
This function creates time series plots of foreign reduction measurements from FLORIDyn simulations. It handles:
1. Time normalization by subtracting the start time
2. Foreign reduction plots in either separated (subplot) or combined layout
3. Different measurement types based on the `msr` parameter

# Plotting Modes
- **Separated mode** (`separated=true`): Creates individual subplots for each turbine
- **Combined mode** (`separated=false`): Plots all turbines on a single figure with different colors

# Example
```julia
using ControlPlots

# Plot velocity reduction for all turbines in combined mode
plotMeasurements(plt, wind_farm, measurements_df, vis; msr=VelReduction)

# Plot added turbulence with separated subplots
plotMeasurements(plt, wind_farm, measurements_df, vis; separated=true, msr=AddedTurbulence)
```

# See Also
- [`plotFlowField`](@ref): For flow field visualization
- [`getMeasurements`](@ref): For generating measurement data
"""
function plotMeasurements(plt, wf::WindFarm, md::DataFrame, vis::Vis; separated=false, msr=VelReduction, pltctrl=nothing)
    # local fig
    # Master switch: still create figures for saving, but suppress interactive pauses/display when show_plots=false
    show_plots = vis.show_plots
    # fig = nothing  # Initialize fig variable

    # Subtract start time
    timeFDyn = md.Time .- md.Time[1]

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
    
    plot_size = 1
    measurement_data = md[!, data_column]
    
    # Measurement plotting
    if separated
        # Calculate y-axis limits
        y_lim = [minimum(measurement_data), maximum(measurement_data)]
        y_range = y_lim[2] - y_lim[1]
        y_lim = y_lim .+ [-0.1, 0.1] * max(y_range, 0.5)
        if wf.nT < 10
            lay = get_layout(wf.nT)
            fig = plt.figure(title, figsize=(10*plot_size, 6*plot_size))
            c = plt.get_cmap("inferno")(0.5)  # Single color for separated plots
            
            for iT in 1:wf.nT
                plt.subplot(lay[1], lay[2], iT)
                plt.plot(
                    timeFDyn[iT:wf.nT:end],
                    measurement_data[iT:wf.nT:end],
                    linewidth=2, color=c
                )
                plt.grid(true)
                plt.title("Turbine $(iT)")
                plt.xlim(max(timeFDyn[1], 0), timeFDyn[end])
                plt.ylim(y_lim...)
                plt.xlabel("Time [s]")
                plt.ylabel(ylabel)
                plt.tight_layout()   
                fig.subplots_adjust(wspace=0.55)
            end
        else          
            rows, lines = get_layout(wf.nT)
            n_turbines = wf.nT

            # Group turbines into subplots based on layout
            times = Float64[]
            measurements = Vector{Float64}[]
            turbines = 1:wf.nT

            # Extract measurement data for each turbine
            measurement_data = md[!, data_column]
            timeFDyn = md.Time .- md.Time[1]

            # Use the actual time data from the simulation results
            times = timeFDyn[1:wf.nT:end]  # Extract times corresponding to first turbine data points

            # Arrays to store time series data  
            measurements = Vector{Float64}[]

            for iT in 1:wf.nT
                push!(measurements, measurement_data[iT:wf.nT:end])
            end

            # Convert vector of vectors to matrix for easier plotting
            # measurements is a vector of 9 vectors, each with 301 time points
            # We want a matrix that's 301 × 9 (time × turbines)
            msr_matrix = hcat(measurements...)  # This creates a 301 × 54 matrix

            # Create dynamic plot arguments based on number of turbines
            n_turbines = wf.nT
            rows, lines = get_layout(wf.nT)

            # Group turbines into subplots based on layout
            plot_data = []
            turbine_labels = []
            subplot_labels = []
                
            turbine_idx = 1
            for row in 1:rows
                # global turbine_idx, lines_in_subplot, labels_in_subplot
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
                
                push!(turbine_labels, ylabel)  # Use the appropriate ylabel for the measurement type
                push!(subplot_labels, labels_in_subplot)
            end
            # Plot with multiple lines per subplot
            try
                plot_x(times, plot_data...; ylabels=turbine_labels, labels=subplot_labels,
                        fig=title, xlabel="rel_time [s]", ysize = 9, bottom=0.02, pltctrl, legend_size=6, loc="center left")
            catch e
                println("Error in plot_x: $e")
            end
        end
    else # not separated
        fig = plt.figure(title*" - Line Plot")
        # Generate colors for each turbine
        color_values = range(0, 1, length=wf.nT)
        colors = [plt.get_cmap("inferno")(val) for val in color_values]
        
        for iT in 1:wf.nT
            plt.plot(
                timeFDyn[iT:wf.nT:end],
                measurement_data[iT:wf.nT:end],
                linewidth=2, color=colors[iT]
            )
        end
        plt.grid(true)
        plt.title(ylabel)
        plt.xlim(max(timeFDyn[1], 0), timeFDyn[end])
        plt.xlabel("Time [s]")
        plt.ylabel(ylabel)
    end
    # plt.show(fig)
    # Save plot to video or output folder if requested
    if vis.save && !vis.unit_test
        directory = vis.output_path
        
        # Generate filename with measurement type and time information
        if separated
            filename = joinpath(directory, "$(msr_name).png")
        else
            filename = joinpath(directory, "$(msr_name)-lineplot.png")
        end
        
        # Save the current figure
        try
            plt.savefig(filename, dpi=150, bbox_inches="tight", pad_inches=0.1, facecolor="white")
            if !vis.unit_test && vis.print_filenames
                @info "Saving $filename"
            end
        catch e
            @warn "Failed to save plot: $e"
        end
    end

    if vis.unit_test
        if show_plots
            plt.pause(1.0)
        end
        close_all(plt)
    end
    if show_plots
        plt.pause(0.01)
    end
    return nothing
end

"""
    get_layout(nT::Int) -> (Int, Int)

Calculate optimal subplot layout (rows, columns) for a given number of turbines.

# Arguments
- `nT::Int`: Number of turbines

# Returns
- `Tuple{Int, Int}`: (rows, columns) for subplot arrangement

# Description
Determines the most square-like arrangement of subplots to accommodate `nT` plots.
"""
function get_layout(nT::Int)
    if nT <= 0
        return (1, 1)
    elseif nT == 1
        return (1, 1)
    elseif nT <= 4
        return (2, 2)
    elseif nT <= 6
        return (2, 3)
    elseif nT <= 9
        return (3, 3)
    elseif nT <= 12
        return (3, 4)
    elseif nT <= 16
        return (4, 4)
    else
        # For larger numbers, calculate a roughly square layout
        cols = ceil(Int, sqrt(nT))
        rows = ceil(Int, nT / cols)
        return (rows, cols)
    end
end

