# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Main script to run a model predictive control (MPC) simulation with FLORIDyn.jl
# Currently, two modes are supported: control of all turbines with the same induction factor using 3 parameter
# and group control with 10 parameters (3 for induction scaling and 7 for individual group scaling).
# The mean square error between the production and demand is minimized.

using Pkg
if ! ("NOMAD" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using FLORIDyn, TerminalPager, DistributedNext, DataFrames, NOMAD, JLD2, Statistics, Printf
using FLORIDyn: TurbineGroup, TurbineArray
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"
data_file               = "data/mpc_result.jld2"
error_file              = "data/mpc_error.jld2"
data_file_group_control = "data/mpc_result_group_control.jld2"

GROUPS = 8
GROUP_CONTROL = true  # if false, use 3-parameter control for all turbines; if true, use 10-parameter group control
SIMULATE = true       # if false, load cached results if available
MAX_STEPS = 1      # maximum number black-box evaluations for NOMAD optimizer
USE_TGC = false
USE_STEP = false
USE_FEED_FORWARD = true # if false, use constant induction (no feed-forward)
ONLINE = false
T_SKIP = 400    # skip first 400s of simulation for error calculation and plotting
T_START = 240   # time to start increasing demand
T_END   = 960   # time to reach final demand
T_EXTRA = 1520  # extra time in addition to sim.end_time for MPC simulation
MAX_DISTANCES = Float64[]
DELTA_P = Float64[]
if isfile(error_file)
    DELTA_P = JLD2.load(error_file)["delta_p"]
end

"""
    create_8_groups(ta::TurbineArray) -> Vector{Dict}

Create 8 turbine groups by dividing turbines based on their X coordinates.
Returns a turbine_groups structure compatible with FLORIDyn.

# Arguments
- `ta::TurbineArray`: The turbine array containing position data

# Returns
- `Vector{Dict}`: Vector of group dictionaries with keys "name", "id", and "turbines"
"""
function create_8_groups(ta::TurbineArray)
    n_turbines = size(ta.pos, 1)
    x_coords = ta.pos[:, 1]
    
    # Create array of (turbine_id, x_coord) pairs
    turbines_with_x = [(i, x_coords[i]) for i in 1:n_turbines]
    
    # Sort by X coordinate
    sort!(turbines_with_x, by = x -> x[2])
    
    # Split into 8 groups
    n_groups = 8
    turbines_per_group = div(n_turbines, n_groups)
    remainder = n_turbines % n_groups
    
    turbine_groups = []
    start_idx = 1
    
    for group_id in 1:n_groups
        # Distribute remainder turbines to first groups
        group_size = turbines_per_group + (group_id <= remainder ? 1 : 0)
        end_idx = start_idx + group_size - 1
        
        # Extract turbine IDs for this group
        group_turbines = [turbines_with_x[i][1] for i in start_idx:end_idx]
        
        push!(turbine_groups, Dict(
            "name" => "group_$group_id",
            "id" => group_id,
            "turbines" => group_turbines
        ))
        
        start_idx = end_idx + 1
    end
    
    # Add "all" group
    push!(turbine_groups, Dict(
        "name" => "all",
        "id" => 0,
        "turbines" => collect(1:n_turbines)
    ))
    
    return turbine_groups
end

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
include("remote_plotting.jl")
include("calc_induction_matrix.jl")

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn, ta = setup(settings_file)

# Override with 8 groups if GROUPS == 8
if GROUPS == 8
    println("Creating 8 turbine groups based on X coordinates...")
    turbine_groups = create_8_groups(ta)
    # Convert to TurbineGroup objects
    new_groups = [TurbineGroup(g["name"], g["id"], g["turbines"]) for g in turbine_groups]
    ta = TurbineArray(ta.pos, ta.type, ta.init_States, new_groups)
end

sim.end_time += T_EXTRA  # extend simulation time for MPC
con.yaw="Constant"
con.yaw_fixed = 270.0
wind.input_dir="Constant"
wind.dir_fixed = 270.0
induction = calc_induction_per_group(1, 0)
set_induction!(ta, induction)

time_step = sim.time_step  # seconds
t_end = sim.end_time - sim.start_time  # relative end time in seconds
con.induction_data = calc_induction_matrix(ta, con, time_step, t_end)

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
if USE_FEED_FORWARD
    set.induction_mode = Induction_TGC()
else
    set.induction_mode = Induction_Constant()
end
wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# Calculate demand for each time point
time_vector = 0:time_step:t_end
demand_values = [calc_demand(t) for t in time_vector]

"""
    calc_max_power(wind_speed, ta, wf, floris) -> Float64

Calculate the theoretical maximum power output for the wind farm.

# Arguments
- `wind_speed::Float64`: Free flow wind speed in m/s
- `ta::TurbineArray`: Turbine array containing position data
- `wf::WindFarm`: Wind farm object containing rotor diameter data
- `floris::Floris`: FLORIS parameters containing air density and efficiency

# Returns
- `max_power::Float64`: Maximum total power in MW

# Assumptions
- Optimal axial induction factor (Betz limit: a = 1/3)
- No yaw misalignment (yaw = 0°)
- All turbines have the same rotor diameter
"""
function calc_max_power(wind_speed, ta, wf, floris)
    a_opt = 1/3  # optimal axial induction factor (Betz limit)
    Cp_opt = 4 * a_opt * (1 - a_opt)^2  # optimal power coefficient
    yaw = 0.0  # no yaw angle

    # Calculate maximum power for all turbines using getPower formula
    # P = 0.5 * ρ * A * Cp * U^3 * η * cos(yaw)^p_p
    nT = length(ta.pos[:, 1])  # number of turbines
    rotor_area = π * (wf.D[1] / 2)^2  # assuming all turbines have same diameter
    max_power_per_turbine = 0.5 * floris.airDen * rotor_area * Cp_opt * wind_speed^3 * floris.eta * cos(yaw)^floris.p_p / 1e6  # MW
    max_power = nT * max_power_per_turbine  # total maximum power in MW
end

# This function implements the "model" in the block diagram.
function run_simulation(set_induction::AbstractMatrix)
    global set, wind, con, floridyn, floris, sim, ta, vis 
    con.induction_data = set_induction
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    vis.online = ONLINE
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    # Calculate total wind farm power by grouping by time and summing turbine powers
    total_power_df = combine(groupby(md, :Time), :PowerGen => sum => :TotalPower)
    # Calculate theoretical maximum power based on turbine ratings and wind conditions
    # Assumptions: free flow wind speed from wind.vel, optimal axial induction factor, no yaw
    max_power = calc_max_power(wind.vel, ta, wf, floris)
    rel_power = (total_power_df.TotalPower ./ max_power) 
end

function calc_axial_induction2(time, scaling::Vector; dt=T_SKIP, group_id=nothing)
    distance = 0.0
    id_scaling = 1.0
    if length(scaling) > 3 && !isnothing(group_id)
        if group_id >= 1 && group_id <= 7
            id_scaling = scaling[3 + group_id]
        elseif group_id == 8
            # Group 8: calculate as 8.0 minus sum of groups 1-7
            id_scaling = 8.0 - sum(scaling[4:10])
        end
        id_scaling = clamp(id_scaling, 0.0, 2.0)
    end
    t1 = 240.0 + dt  # Time to start increasing demand
    t2 = 960.0 + dt  # Time to reach final demand

    if time < t1
        time = t1
    end

    scaling_begin = scaling[1]
    scaling_mid = scaling[2]
    scaling_end = scaling[3]
    
    # Monotonic cubic interpolation using Fritsch-Carlson method
    # This prevents overshoot/undershoot while maintaining smoothness
    s = clamp((time - t1) / (t2 - t1), 0.0, 1.0)
    
    # Three control points at s = 0, 0.5, 1
    # Using piecewise quadratic Bezier curves for true monotonicity
    if s <= 0.5
        # First segment: quadratic interpolation from begin to mid
        t_local = s / 0.5  # normalize to [0, 1]
        # Quadratic Bezier: P(t) = (1-t)^2*P0 + 2t(1-t)*P1 + t^2*P2
        # Choose P1 as linear interpolation to ensure monotonicity
        p0 = scaling_begin
        p2 = scaling_mid
        p1 = 0.5 * (p0 + p2)  # midpoint ensures no overshoot
        scaling_result = (1 - t_local)^2 * p0 + 2 * t_local * (1 - t_local) * p1 + t_local^2 * p2
    else
        # Second segment: quadratic interpolation from mid to end
        t_local = (s - 0.5) / 0.5  # normalize to [0, 1]
        p0 = scaling_mid
        p2 = scaling_end
        p1 = 0.5 * (p0 + p2)  # midpoint ensures no overshoot
        scaling_result = (1 - t_local)^2 * p0 + 2 * t_local * (1 - t_local) * p1 + t_local^2 * p2
    end
    
    demand = calc_demand(time)
    demand_end = calc_demand(t2)
    interpolated_demand = demand_end - (demand_end - demand) * id_scaling
    scaled_demand = scaling_result * interpolated_demand
    if scaled_demand > 1.0
        distance = scaled_demand - 1.0
        # @warn("Scaled demand exceeds 100% at time=$(time)s: scaled_demand=$(scaled_demand), scaling_result=$(scaling_result), demand=$(demand), id_scaling=$(id_scaling)")
    end
    base_induction = calc_induction(scaled_demand * cp_max)

    # Calculate interpolation factor
    # 1.0 at t=t1 (full correction), 0.0 at t=t2 (no correction)
    # Linear interpolation between t1 and t2
    if time <= t1
        interp_factor = 1.0  # Full correction before and at t1
    elseif time >= t2
        interp_factor = 0.0  # No correction at and after t2
    else
        # Linear interpolation between t1 and t2
        interp_factor = ((t2 - time) / (t2 - t1))
    end
    
    rel_power = calc_cp(base_induction) / cp_max
    if length(DELTA_P) == 0
        delta_p = 0.0
    else
        delta_p = DELTA_P[Int(floor((time) / time_step)) + 1]  # +1 for 1-based indexing
    end
    if isnan(delta_p)
        delta_p = 0.0
    end
    corrected_rel_power = rel_power - 0.0000000000001 * delta_p
    corrected_induction = calc_induction(corrected_rel_power * cp_max)
    return max(0.0, min(BETZ_INDUCTION, corrected_induction)), distance
end

function calc_induction_matrix2(ta, time_step, t_end; scaling)
    # Create time vector from 0 to t_end with time_step intervals
    time_vector = 0:time_step:t_end
    n_time_steps = length(time_vector)
    n_turbines = size(ta.pos, 1)  # Use ta.pos to get number of turbines
    
    # Initialize matrix: rows = time steps, columns = time + turbines
    # First column is time, subsequent columns are turbine induction factors
    induction_matrix = zeros(Float64, n_time_steps, n_turbines + 1)
    
    # Fill the first column with time values
    induction_matrix[:, 1] = collect(time_vector)
    
    # Calculate induction for each turbine at each time step (columns 2 onwards)
    max_distance = 0.0
    for (t_idx, time) in enumerate(time_vector)
        for i in 1:n_turbines
            group_id = FLORIDyn.turbine_group(ta, i)
            axial_induction, distance = calc_axial_induction2(time, scaling; group_id=group_id)
            induction_matrix[t_idx, i + 1] = axial_induction
            max_distance = max(max_distance, distance)
        end
    end

    return induction_matrix, max_distance
end

function calc_error(rel_power, demand_values, time_step)
    # Start index after skipping initial transient; +1 because Julia is 1-based
    i0 = Int(floor(T_SKIP / time_step)) + 1
    # Clamp to valid range
    i0 = max(1, i0)
    n = min(length(rel_power), length(demand_values)) - i0 + 1
    if n <= 0
        error("calc_error: empty overlap after skip; check T_SKIP and lengths (rel_power=$(length(rel_power)), demand=$(length(demand_values)), i0=$(i0))")
    end
    r = @view rel_power[i0:i0 + n - 1]
    d = @view demand_values[i0:i0 + n - 1]
    return sum((r .- d) .^ 2) / length(d)
end

"""
    plot_induction(optimal_scaling::Vector{Float64})

Plot the axial induction factor for turbine group 1 over time range 500-1500s.

# Arguments
- `optimal_scaling::Vector{Float64}`: Optimal scaling parameters from optimization

# Description
Uses `calc_axial_induction2` to compute induction values for group 1 and plots
them against time. The plot uses the global `pltctrl` variable for thread-safe plotting.
"""
function plot_induction(optimal_scaling::Vector{Float64})
    # Time range: 500 to 1500 seconds
    t_start = 500.0
    t_end = 1500.0
    dt = time_step
    
    # Print diagnostic information
    println("\n=== Diagnostic: plot_induction ===")
    println("optimal_scaling[1:3]: ", optimal_scaling[1:3])
    if length(optimal_scaling) > 3
        println("optimal_scaling[4] (group 1 id_scaling): ", optimal_scaling[4])
    end
    
    # Create time vector
    time_vec = t_start:dt:t_end
    n_points = length(time_vec)
    
    # Calculate induction for group 1 at each time point
    induction_values = zeros(n_points)
    group_id = 1
    
    for (i, t) in enumerate(time_vec)
        induction, _ = calc_axial_induction2(t, optimal_scaling; group_id=group_id)
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

# Calculate storage time at 100% power in seconds
function calc_storage_time(time_vector, rel_power_gain)
    dt = time_vector[2]-time_vector[1]
    mean_gain = mean(rel_power_gain)
    time = length(rel_power_gain)*dt
    storage_time = mean_gain * time
end

"""
    eval_fct(x::Vector{Float64}) -> Tuple{Bool, Bool, Vector{Float64}}

Evaluation function for NOMAD optimization of the scaling parameter.

This function calculates the mean squared error between the relative power output 
and demand values for a given scaling parameter. It is designed to be minimized 
by the NOMAD optimizer.

# Arguments
- `x::Vector{Float64}`: A vector containing the scaling parameters

# Returns
- `success::Bool`: Always `true` to indicate successful evaluation
- `count_eval::Bool`: Always `true` to count this evaluation
- `bb_outputs::Vector{Float64}`: Vector containing [objective, constraint(s)]

# Global Variables Used
- `ta`: Turbine array
- `time_step`: Simulation time step
- `t_end`: End time of simulation
- `demand_values`: Target demand values
- `GROUP_CONTROL`: Boolean flag for group control mode
"""
function eval_fct(x::Vector{Float64})
    scaling = x  # scaling is now a vector with parameters
    print(".")  # progress indicator
    
    # Calculate induction matrix with current scaling
    induction_data, max_distance = calc_induction_matrix2(ta, time_step, t_end; scaling=scaling)
    if max_distance > 0.0
        push!(MAX_DISTANCES, max_distance)
    end

    # Run simulation and get relative power
    rel_power = run_simulation(induction_data)
    
    # Calculate error
    error = calc_error(rel_power, demand_values, time_step)
    
    # Add constraint if GROUP_CONTROL is true
    if GROUP_CONTROL
        # Constraint: x[4] + x[5] + ... + x[10] <= 8
        # For NOMAD, constraints should be <= 0, so we formulate as:
        # x[4] + x[5] + ... + x[10] - 8 <= 0
        constraint_sum = sum(x[4:10]) - 8.0
        # Constraint 2: max_distance <= 0.075
        # Formulate as: max_distance - 0.075 <= 0
        constraint_maxdist = max_distance - 0.075*3
        bb_outputs = [error, constraint_sum, constraint_maxdist]
    else
        bb_outputs = [error]
    end
    
    success = true
    count_eval = true
    
    return (success, count_eval, bb_outputs)
end
if GROUP_CONTROL
    # Set up NOMAD optimization problem
    p = NomadProblem(
        10,                   # dimension (10 parameters: scaling_begin, scaling_mid, scaling_end, id_scaling for groups 1-7)
        3,                    # number of outputs (objective + 2 constraints)
        ["OBJ", "PB", "PB"], # output types: OBJ = objective to minimize, PB = progressive barrier constraints
        eval_fct;            # evaluation function
        lower_bound=[1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # minimum scaling values (3 global + 7 groups)
        upper_bound=[2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]    # maximum scaling values
    )

    # Set NOMAD options
    p.options.max_bb_eval = MAX_STEPS      # maximum number of function evaluations
    p.options.display_degree = 2    # verbosity level
else
        # Set up NOMAD optimization problem
    p = NomadProblem(
        3,                    # dimension (3 parameters: scaling_begin, scaling_mid, scaling_end)
        1,                    # number of outputs (just the objective)
        ["OBJ"],             # output types: OBJ = objective to minimize
        eval_fct;            # evaluation function
        lower_bound=[1.0, 1.0, 1.0],   # minimum scaling values
        upper_bound=[2.0, 3.0, 3.0]    # maximum scaling values
    )

    # Set NOMAD options
    p.options.max_bb_eval = MAX_STEPS      # maximum number of function evaluations
    p.options.display_degree = 2           # verbosity level
end

results = nothing
if (! SIMULATE) && ((isfile(data_file) && !GROUP_CONTROL) || (isfile(data_file_group_control) && GROUP_CONTROL))
    println("Loading cached MPC results from $(data_file)…")
    if GROUP_CONTROL
        results = JLD2.load(data_file_group_control, "results")
        results_ref = JLD2.load(data_file, "results")
        rel_power_ref = results_ref["rel_power"]
    else
        results = JLD2.load(data_file, "results")   
    end
    # Unpack
    time_vector   = results["time_vector"]
    demand_values = results["demand_values"]
    rel_power     = results["rel_power"]
    induction_data = results["induction_data"]
    optimal_scaling = results["optimal_scaling"]
    mse = results["mse"]
else
    # Run optimization and simulation
    if GROUP_CONTROL
        #  [1.25773, 1.25069, 1.27628, 0.000114, 0.00029, 1.99977, 1.99912, 1.5639, 0.00031, 0.4386]
        result = solve(p, [1.25773, 1.25069, 1.27628, 0.000114, 0.00029, 1.99977, 1.99912, 1.5639, 0.00031, 0.4386])
        results_ref = JLD2.load(data_file, "results")
        rel_power_ref = results_ref["rel_power"]
        optimal_scaling = result.x_best_feas[1:10]
    else
        result = solve(p, [1.5, 1.5, 1.5])  # Start from initial guess of [1.5, 1.5, 1.5]
        optimal_scaling = result.x_best_feas[1:3]
    end

    induction_data, max_distance = calc_induction_matrix2(ta, time_step, t_end; scaling=optimal_scaling)
    rel_power = run_simulation(induction_data)
    mse = calc_error(rel_power, demand_values, time_step)

    # Persist
    if GROUP_CONTROL
        data_file1 = data_file_group_control
    else
        data_file1 = data_file
    end
    JLD2.jldsave(data_file1; results=Dict(
        "time_vector" => collect(time_vector),
        "demand_values" => demand_values,
        "rel_power" => rel_power,
        "induction_data" => induction_data,
        "optimal_scaling" => optimal_scaling,
        "mse" => mse,
    ))
end

println("\nRoot Mean Square Error (RMSE): $(round(sqrt(mse) * 100, digits=2))%")
# calculate delta_p
delta_p = rel_power[1:length(time_vector)].-demand_values
if length(DELTA_P) == 0
    JLD2.jldsave(error_file; delta_p=delta_p, time_vector=time_vector)
end

if GROUP_CONTROL
    plot_rmt(time_vector, [rel_power[1:length(time_vector)] .* 100, rel_power_ref[1:length(time_vector)] .* 100, demand_values .* 100]; xlabel="Time [s]", xlims=(T_SKIP, time_vector[end]),
            ylabel="Rel. Power Output [%]", labels=["rel_power", "rel_power_ref", "rel_demand"], fig="Rel. Power and Demand", pltctrl)
else
    plot_rmt(time_vector, [rel_power[1:length(time_vector)] .* 100, demand_values .* 100]; xlabel="Time [s]", xlims=(T_SKIP, time_vector[end]),
            ylabel="Rel. Power Output [%]", labels=["rel_power", "rel_demand"], fig="Rel. Power and Demand", pltctrl)
end
# ## plot induction factor vs time for one turbine using calc_axial_induction2
# induction_factors = induction_data[:, 2]
# plot_rmt(time_vector, induction_factors; xlabel="Time [s]", ylabel="Axial Induction Factor", fig="induction", pltctrl)

# Plot average axial induction factor per turbine group over time
begin
    # Extract time from first column of induction_data
    time_vec_ind = induction_data[:, 1]
    n_time_steps = size(induction_data, 1)
    n_turbines = size(ta.pos, 1)

    # Prepare containers for eight groups
    group_data = [Float64[] for _ in 1:8]

    # Since all turbines in a group have identical induction, pick one representative turbine per group
    group_indices = [findfirst(i -> FLORIDyn.turbine_group(ta, i) == g, 1:n_turbines) for g in 1:8]
    for g in 1:8
        idx = group_indices[g]
        if isnothing(idx)
            group_data[g] = fill(0.0, n_time_steps)
        else
            group_data[g] = induction_data[:, idx + 1]
        end
    end

    group_labels = ["Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8"]
    plot_rmt(time_vec_ind, group_data;
             xlabel="Time [s]",
             ylabel="Axial Induction Factor [-]",
             title="Average Axial Induction Factor vs Time by Turbine Group",
             labels=group_labels,
             fig="Induction by Group",
             pltctrl=pltctrl)
    
    # Create bar plot of average induction factor per group
    avg_induction = [isempty(group_data[g]) ? 0.0 : mean(group_data[g]) for g in 1:8]
    
    if !isnothing(plt)
        plt.figure(figsize=(10, 6))
        bars = plt.bar(1:8, avg_induction, color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"])
        plt.xlabel("Turbine Group", fontsize=12)
        plt.ylabel("Average Axial Induction Factor [-]", fontsize=12)
        plt.title("Average Axial Induction Factor by Turbine Group", fontsize=14)
        plt.xticks(1:8, ["Group $i" for i in 1:8], rotation=45, ha="right")
        plt.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        
        # Add value labels on top of bars
        for (i, v) in enumerate(avg_induction)
            plt.text(i, v, @sprintf("%.3f", v), ha="center", va="bottom", fontsize=10)
        end
    end
end

if GROUP_CONTROL
    # calculate rel_power-rel_power_ref
    start_index = Int(floor((T_SKIP+T_START+(T_END-T_START)*0.96) / time_step)) + 1
    rel_power_gain = rel_power[start_index:end] .- rel_power_ref[start_index:end]
    storage_time = calc_storage_time(time_vector, rel_power_gain)
    println("Estimated storage time at 100% power: $(round(storage_time, digits=2)) s")
    println()
    plot_rmt((1:length(rel_power_gain)).*4, rel_power_gain .* 100; xlabel="Time [s]", ylabel="Rel. Power Gain [%]", fig="rel_power_ref", pltctrl)
    results = JLD2.load(data_file_group_control, "results")
else
    results = JLD2.load(data_file, "results")
end

