# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Main script to run a model predictive control (MPC) simulation with FLORIDyn.jl
# Currently, two modes are supported: control of all turbines with the same induction factor using CONTROL_POINTS parameters
# and group control with CONTROL_POINTS + (GROUPS-1) parameters (CONTROL_POINTS for induction scaling at different time points 
# and GROUPS-1 for individual group scaling, with the last group calculated from a constraint).

# The script uses the NOMAD.jl package for black-box optimization of the scaling parameters.
# The number of groups can be 1, 2, 3, 4, 8, or 12.
# The mean square error between the production and demand is minimized.

# Always run the script with GROUPS = 1 first to get a baseline result without group control. This baseline 
# is stored in the file data/mpc_result.jld2 and used for comparison when running with group control.

# The constant MAX_STEPS can be one for debugging. With 10 or 100 you already get a rough idea of the optimization 
# progress. So far, I never needed more than about 1600 steps to converge with GROUPS = 12 and less with fewer groups.

# If you want to create a video, make sure to set ONLINE = true and run Julia with sufficient threads. On a 
# 7850X with 16 threads, the video creation needs about two hours.

# To create a bar plot, run Julia single threaded.

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
data_file_group_control = "data/mpc_result_group_control"

GROUPS = 3 # must be 1, 2, 3, 4, 8 or 12
CONTROL_POINTS = 5
MAX_ID_SCALING = 3.0
SIMULATE = true      # if false, load cached results if available
MAX_STEPS = 4000       # maximum number black-box evaluations for NOMAD optimizer
USE_TGC = false
USE_STEP = false
USE_FEED_FORWARD = true # if false, use constant induction (no feed-forward)
ONLINE  = false  # if true, enable online plotting during simulation and create video
T_START = 240   # relative time to start increasing demand
T_END   = 960   # relative time to reach final demand
T_EXTRA = 2580  # extra time in addition to sim.end_time for MPC simulation
MIN_INDUCTION = 0.01
MAX_DISTANCES = Float64[]
data_file_group_control = data_file_group_control * '_' * string(GROUPS) * "TGs.jld2"

GROUP_CONTROL = (GROUPS != 1)
@assert(GROUPS in (1, 2, 3, 4, 8, 12), "GROUPS must be 1, 2, 3, 4, 8, or 12")

# Load vis settings from YAML file
vis = Vis(vis_file)
vis.save = ONLINE
# For GROUP_CONTROL, disable online visualization during initial setup to avoid NaN issues
# The visualization will be enabled after the first valid induction_data is calculated
vis.online = ONLINE && !GROUP_CONTROL
if ONLINE
    cleanup_video_folder()
end
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

# Override with n groups if GROUPS != 4 (default in settings file is 4)
if GROUPS != 4
    println("Creating $GROUPS turbine groups based on X coordinates...")
    ta = create_n_groups(ta, GROUPS)
end

sim.end_time += T_EXTRA  # extend simulation time for MPC
con.yaw="Constant"
con.yaw_fixed = 270.0
wind.input_dir="Constant"
wind.dir_fixed = 270.0
induction = calc_induction_per_group(vis, 1, 0)
set_induction!(ta, induction)

time_step = sim.time_step  # seconds
t_end = sim.end_time - sim.start_time  # relative end time in seconds

# For initial setup, use calc_induction_matrix (only affects pre-optimization visualization)
# During optimization, calc_induction_matrix2 will be used with proper group handling
con.induction_data = calc_induction_matrix(vis, ta, time_step, t_end)

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
demand_values = [calc_demand(vis, t) for t in time_vector]

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
function run_simulation(set_induction::AbstractMatrix; enable_online=false)
    global set, wind, con, floridyn, floris, sim, ta, vis 
    con.induction_data = set_induction
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    # Only enable online visualization if explicitly requested (to avoid NaN issues during optimization)
    vis.online = enable_online
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
    # Calculate total wind farm power by grouping by time and summing turbine powers
    total_power_df = combine(groupby(md, :Time), :PowerGen => sum => :TotalPower)
    # Calculate theoretical maximum power based on turbine ratings and wind conditions
    # Assumptions: free flow wind speed from wind.vel, optimal axial induction factor, no yaw
    max_power = calc_max_power(wind.vel, ta, wf, floris)
    rel_power = (total_power_df.TotalPower ./ max_power) 
end

"""
    interpolate_hermite_spline(s::Float64, scaling::Vector) -> Float64

Perform piecewise cubic Hermite spline interpolation across control points.

This function uses cubic Hermite spline interpolation between control points 
evenly spaced along s ∈ [0, 1]. The method provides smooth C1-continuous transitions while
respecting the control point values. The number of control points is determined
from the length of the `scaling` vector.

# Arguments
- `s::Float64`: Normalized parameter in [0, 1] representing position along the curve
- `scaling::Vector`: Vector containing control point values. For n control points,
  they are located at s = 0, 1/(n-1), 2/(n-1), ..., 1.0

# Returns
- `Float64`: Interpolated value at position s
"""
function interpolate_hermite_spline(s::Float64, scaling::Vector)
    n_points = length(scaling)
    @assert n_points >= 2 "Need at least 2 control points for interpolation"
    
    # Handle edge cases
    if n_points == 2
        # Linear interpolation for 2 points
        return scaling[1] + s * (scaling[2] - scaling[1])
    end
    
    # Number of segments = number of control points - 1
    n_segments = n_points - 1
    segment_width = 1.0 / n_segments
    
    # Calculate tangents for all control points using finite differences
    tangents = zeros(n_points)
    
    # First point: forward difference
    tangents[1] = (scaling[2] - scaling[1]) / segment_width
    
    # Interior points: central differences
    for i in 2:(n_points-1)
        tangents[i] = (scaling[i+1] - scaling[i-1]) / (2 * segment_width)
    end
    
    # Last point: backward difference
    tangents[n_points] = (scaling[n_points] - scaling[n_points-1]) / segment_width
    
    # Determine which segment we're in
    segment_idx = min(n_segments, Int(floor(s / segment_width)) + 1)
    if s >= 1.0
        segment_idx = n_segments
    end
    
    # Local parameter t within the segment [0, 1]
    s_start = (segment_idx - 1) * segment_width
    t = (s - s_start) / segment_width
    t = clamp(t, 0.0, 1.0)
    
    # Cubic Hermite basis functions
    h00 = 2*t^3 - 3*t^2 + 1
    h10 = t^3 - 2*t^2 + t
    h01 = -2*t^3 + 3*t^2
    h11 = t^3 - t^2
    
    # Interpolate using values and tangents at segment endpoints
    p0 = scaling[segment_idx]
    p1 = scaling[segment_idx + 1]
    m0 = tangents[segment_idx]
    m1 = tangents[segment_idx + 1]
    
    scaling_result = h00*p0 + h10*segment_width*m0 + h01*p1 + h11*segment_width*m1
    
    return scaling_result
end

"""
    calc_axial_induction2(vis, time, scaling::Vector; group_id=nothing) -> (corrected_induction, distance)

Calculate the axial induction factor for a turbine using optimizable scaling parameters.

This function computes the axial induction factor based on:
1. Time-dependent scaling via cubic Hermite spline interpolation of control points
2. Optional group-specific scaling factors for individual turbine group control
3. Demand-based adjustment with correction for power coefficient nonlinearity

# Arguments
- `vis`: [`Vis`](@ref) object containing visualization settings (uses `t_skip`)
- `time::Float64`: Current simulation time in seconds
- `scaling::Vector{Float64}`: Optimization parameters vector containing:
  - Elements 1 to CONTROL_POINTS: Time-dependent scaling control points
  - Elements CONTROL_POINTS+1 to end: Group-specific scaling factors (if GROUP_CONTROL)
- `group_id::Union{Int,Nothing}`: Turbine group identifier (1 to GROUPS), or `nothing` for no group control

# Returns
- `corrected_induction::Float64`: Computed axial induction factor, clamped to [MIN_INDUCTION, BETZ_INDUCTION]
- `distance::Float64`: Constraint violation distance (positive if scaled_demand > 1.0, else 0.0)

# Details
The function operates in several stages:
1. Extracts group-specific scaling factor `id_scaling` from the scaling vector (if applicable)
   - For groups 1 to GROUPS-1: directly from `scaling[CONTROL_POINTS + group_id]`
   - For the last group (GROUPS): calculated as `GROUPS * MAX_ID_SCALING / 2.0 - sum(scaling[(CONTROL_POINTS+1):end])`
     to reduce the number of optimization variables by one
2. Computes normalized time parameter `s` ∈ [0,1] between T_START and T_END
3. Interpolates time-dependent scaling using [`interpolate_hermite_spline`](@ref)
4. Adjusts demand by group-specific scaling and applies time-dependent scaling
5. Converts scaled demand to induction, applies power coefficient correction
6. Ensures minimum induction to avoid numerical issues in wake model

# Global Constants Used
- `CONTROL_POINTS`: Number of time-dependent control points
- `GROUPS`: Number of turbine groups
- `MAX_ID_SCALING`: Maximum allowed group scaling factor
- `T_START`: Time offset to start ramping demand (relative to `vis.t_skip`)
- `T_END`: Time offset to reach final demand (relative to `vis.t_skip`)
- `MIN_INDUCTION`: Minimum induction to prevent NaN in FLORIS wake model
- `BETZ_INDUCTION`: Maximum theoretical induction (Betz limit)

# See Also
- [`interpolate_hermite_spline`](@ref): Performs cubic Hermite spline interpolation
- [`calc_induction_matrix2`](@ref): Uses this function to build induction matrices
"""
function calc_axial_induction2(vis, time, scaling::Vector; group_id=nothing)
    distance = 0.0
    id_scaling = 1.0
    if length(scaling) > CONTROL_POINTS && !isnothing(group_id) && group_id >= 1
        if group_id <= GROUPS - 1
            id_scaling = scaling[CONTROL_POINTS + group_id]
        elseif group_id == GROUPS
            # Last group: calculate as GROUPS * MAX_ID_SCALING / 2.0 minus sum of groups 1 to GROUPS-1
            id_scaling = GROUPS * MAX_ID_SCALING / 2.0 - sum(scaling[(CONTROL_POINTS+1):end])
        end
        id_scaling = clamp(id_scaling, 0.0, MAX_ID_SCALING)
    end
    t1 = vis.t_skip + T_START  # Time to start increasing demand
    t2 = vis.t_skip + T_END    # Time to reach final demand

    if time < t1
        time = t1
    end
    
    s = clamp((time - t1) / (t2 - t1), 0.0, 1.0)
    
    # Perform piecewise cubic Hermite spline interpolation
    scaling_result = interpolate_hermite_spline(s, scaling[1:CONTROL_POINTS])
    
    demand = calc_demand(vis, time)
    demand_end = calc_demand(vis, t2)
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
    corrected_induction = calc_induction(rel_power * cp_max)
    
    # Ensure minimum induction to avoid numerical issues in FLORIS (NaN from zero induction)
    # Minimum value of MIN_INDUCTION ensures the wake model has valid inputs
    corrected_induction = max(MIN_INDUCTION, min(BETZ_INDUCTION, corrected_induction))
    
    return corrected_induction, distance
end

function calc_induction_matrix2(vis, ta, time_step, t_end; scaling)
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
            axial_induction, distance = calc_axial_induction2(vis, time, scaling; group_id=group_id)
            induction_matrix[t_idx, i + 1] = axial_induction
            max_distance = max(max_distance, distance)
        end
    end

    return induction_matrix, max_distance
end

function calc_error(vis, rel_power, demand_values, time_step)
    # Start index after skipping initial transient; +1 because Julia is 1-based
    i0 = Int(floor(vis.t_skip / time_step)) + 1
    # Clamp to valid range
    i0 = max(1, i0)
    n = min(length(rel_power), length(demand_values)) - i0 + 1
    if n <= 0
        error("calc_error: empty overlap after skip; check vis.t_skip and lengths (rel_power=$(length(rel_power)), demand=$(length(demand_values)), i0=$(i0))")
    end
    r = @view rel_power[i0:i0 + n - 1]
    d = @view demand_values[i0:i0 + n - 1]
    return sum((r .- d) .^ 2) / length(d)
end

"""
    plot_induction(vis, optimal_scaling::Vector{Float64})

Plot the axial induction factor for turbine group 1 over time range 500-1500s.

# Arguments
- `vis`: Visualization object containing `t_skip` parameter
- `optimal_scaling::Vector{Float64}`: Optimal scaling parameters from optimization

# Description
Uses `calc_axial_induction2` to compute induction values for group 1 and plots
them against time. The plot uses the global `pltctrl` variable for thread-safe plotting.
"""
function plot_induction(vis, optimal_scaling::Vector{Float64})
    # Time range: 500 to 1500 seconds
    t_start = vis.t_skip
    t_end   = vis.t_skip + T_END + T_EXTRA
    dt = time_step
    
    # Print diagnostic information
    println("\n=== Diagnostic: plot_induction ===")
    println("optimal_scaling[1:$CONTROL_POINTS]: ", optimal_scaling[1:CONTROL_POINTS])
    if length(optimal_scaling) > CONTROL_POINTS
        println("optimal_scaling[$(CONTROL_POINTS+1)] (group 1 id_scaling): ", optimal_scaling[CONTROL_POINTS+1])
    end
    
    # Create time vector
    time_vec = t_start:dt:t_end
    n_points = length(time_vec)
    
    # Calculate induction for group 1 at each time point
    induction_values = zeros(n_points)
    group_id = 1
    
    for (i, t) in enumerate(time_vec)
        induction, _ = calc_axial_induction2(vis, t, optimal_scaling; group_id=group_id)
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
    plot_scaling_curve(optimal_scaling::Vector{Float64})

Plot the scaling curve from the piecewise cubic Hermite spline interpolation over s=0..1.

# Arguments
- `optimal_scaling::Vector{Float64}`: Optimal scaling parameters from optimization

# Description
Plots the piecewise cubic Hermite spline interpolation curve showing how the scaling
factor varies across the normalized parameter s from 0 to 1. Uses the first CONTROL_POINTS
elements of `optimal_scaling` as control points evenly spaced from s = 0 to s = 1.0.
"""
function plot_scaling_curve(optimal_scaling::Vector{Float64})
    # Create s vector from 0 to 1
    s_vec = 0.0:0.01:1.0
    n_points = length(s_vec)
    
    # Calculate scaling_result for each s value
    scaling_values = zeros(n_points)
    
    for (i, s) in enumerate(s_vec)
        scaling_values[i] = interpolate_hermite_spline(s, optimal_scaling[1:CONTROL_POINTS])
    end
    
    # Print diagnostic information
    println("\n=== Diagnostic: plot_scaling_curve ===")
    println("Control points (scaling[1:$CONTROL_POINTS]): ", optimal_scaling[1:CONTROL_POINTS])
    println("Min scaling: $(round(minimum(scaling_values), digits=4))")
    println("Max scaling: $(round(maximum(scaling_values), digits=4))")
    println("======================================\n")
    
    # Plot
    plot_rmt(collect(s_vec), scaling_values;
             xlabel="Normalized Parameter s [-]",
             ylabel="Scaling Factor [-]",
             title="Hermite Spline Interpolation Scaling Curve",
             fig="Scaling Curve",
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
    induction_data, max_distance = calc_induction_matrix2(vis, ta, time_step, t_end; scaling=scaling)
    if max_distance > 0.0
        push!(MAX_DISTANCES, max_distance)
    end

    # Run simulation and get relative power
    rel_power = run_simulation(induction_data)
    
    # Calculate error
    error = calc_error(vis, rel_power, demand_values, time_step)
    success = true
    if isnothing(error) || isnan(error)
        error = 1e6
        success = false
    end
    
    # Add constraint if GROUP_CONTROL is true
    if GROUP_CONTROL
        # Constraint: x[CONTROL_POINTS+1] + x[CONTROL_POINTS+2] + ... <= GROUPS * MAX_ID_SCALING / 2.0
        # For NOMAD, constraints should be <= 0, so we formulate as:
        # x[CONTROL_POINTS+1] + x[CONTROL_POINTS+2] + ... - GROUPS * MAX_ID_SCALING / 2.0 <= 0
        constraint_sum = sum(x[(CONTROL_POINTS+1):end]) - (GROUPS * MAX_ID_SCALING / 2.0)
        bb_outputs = [error, constraint_sum]
    else
        bb_outputs = [error]
    end
    
    count_eval = true
    
    return (success, count_eval, bb_outputs)
end
if GROUP_CONTROL
    n_group_params = GROUPS - 1  # One less because last group is calculated from constraint
    n_total_params = CONTROL_POINTS + n_group_params  # CONTROL_POINTS global scaling + (GROUPS-1) group scaling
    
    # Create lower and upper bounds dynamically
    lower_bound = vcat(fill(1.0, CONTROL_POINTS), fill(0.0, n_group_params))
    upper_bound = vcat(fill(2.0, CONTROL_POINTS), fill(MAX_ID_SCALING, n_group_params))
    
    # Set up NOMAD optimization problem
    p = NomadProblem(
        n_total_params,      # dimension (CONTROL_POINTS global + GROUPS-1 group parameters)
        2,                   # number of outputs (objective + 1 constraint)
        ["OBJ", "PB"],       # output types: OBJ = objective to minimize, PB = progressive barrier constraint
        eval_fct;            # evaluation function
        lower_bound=lower_bound,
        upper_bound=upper_bound
    )

    # Set NOMAD options
    p.options.max_bb_eval = MAX_STEPS      # maximum number of function evaluations
    p.options.display_degree = 2    # verbosity level
else
        # Set up NOMAD optimization problem
    p = NomadProblem(
        CONTROL_POINTS,      # dimension (CONTROL_POINTS parameters: scaling at CONTROL_POINTS time points)
        1,                   # number of outputs (just the objective)
        ["OBJ"],             # output types: OBJ = objective to minimize
        eval_fct;            # evaluation function
        lower_bound=fill(1.0, CONTROL_POINTS),   # minimum scaling values
        upper_bound=vcat([2.5], fill(3.0, CONTROL_POINTS - 1))    # maximum scaling values
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
        # Create initial guess: CONTROL_POINTS global parameters + (GROUPS-1) group parameters
        if GROUPS == 8
            x0 = [1.31, 1.4427, 1.35654, 1.28725, 1.28105, 0.0027, 0.0294, 1.8695, 2.0157, 1.8563, 1.1908, 0.0825]
        elseif GROUPS == 4
            x0 = [1.578, 1.991, 1.54259, 1.33791, 1.27339, 0.017865, 0.886214, 2.87895]
        elseif GROUPS == 2
            x0 = [1.52628, 1.9693, 1.4923, 1.35422, 1.26623, 0.5599]
        elseif GROUPS == 3
            x0 = [1.9, 2.0, 1.7, 1.399, 1.3, 0.05, 1.48]
        elseif GROUPS == 12
            # CONTROL_POINTS global + 11 group parameters (last group calculated from constraint)
            x0 = [1.409, 1.60396, 1.43527, 1.30722, 1.26675, 0.0877, 0.1621, 0.1235, 1.99722, 0.016, 1.9725, 1.34014, 1.8945, 0.85491, 2.8402, 2.0101]
        else
            # Generic initial guess for other group counts
            x0 = vcat(fill(1.5, CONTROL_POINTS), fill(1.0, GROUPS - 1))
        end
        result = solve(p, x0)
        results_ref = JLD2.load(data_file, "results")
        rel_power_ref = results_ref["rel_power"]
        optimal_scaling = result.x_best_feas
    else
        result = solve(p, [1.18291, 1.19575, 1.21248, 1.2409, 1.30345])  # Start from initial guess
        optimal_scaling = result.x_best_feas[1:CONTROL_POINTS]
    end

    induction_data, max_distance = calc_induction_matrix2(vis, ta, time_step, t_end; scaling=optimal_scaling)
    
    # Enable online visualization for the final simulation with optimized parameters
    enable_viz = ONLINE && GROUP_CONTROL
    
    rel_power = run_simulation(induction_data; enable_online=enable_viz)
    mse = calc_error(vis, rel_power, demand_values, time_step)

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

if GROUP_CONTROL
    plot_rmt(time_vector, [rel_power[1:length(time_vector)] .* 100, rel_power_ref[1:length(time_vector)] .* 100, demand_values .* 100]; xlabel="Time [s]", xlims=(vis.t_skip, time_vector[end]),
            ylabel="Rel. Power Output [%]", labels=["rel_power", "rel_power_ref", "rel_demand"], title="Rel. Power and Demand "*string(GROUPS)*" TGs", fig="Rel. Power and Demand", pltctrl)
else
    plot_rmt(time_vector, [rel_power[1:length(time_vector)] .* 100, demand_values .* 100]; xlabel="Time [s]", xlims=(vis.t_skip, time_vector[end]),
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

function print_gains(optimal_scaling)
    if !GROUP_CONTROL || GROUPS == 1
        println("\n=== Power Gain per Turbine Group ===")
        println("Group gains not applicable (GROUP_CONTROL is false or only one group).")
        return
    end
    scaling = optimal_scaling[(CONTROL_POINTS+1):end]
    id_scaling = GROUPS * MAX_ID_SCALING / 2.0 - sum(scaling)
    push!(scaling, id_scaling)
    println("\n=== Power Gain per Turbine Group ===")
    for (i, gain) in enumerate(scaling)
        println("Group $i: $(round(gain, digits=2))")
    end
    println("mean: $(round(mean(scaling), digits=2))")
end

if GROUP_CONTROL
    # calculate rel_power-rel_power_ref
    start_index = Int(floor((vis.t_skip-40+T_START+(T_END-T_START)) / time_step)) + 1
    common_length = min(length(rel_power), length(rel_power_ref))
    rel_power = rel_power[1:common_length]
    rel_power_ref = rel_power_ref[1:common_length]
    rel_power_gain = rel_power[start_index:end] .- rel_power_ref[start_index:end]
    storage_time = calc_storage_time(time_vector, rel_power_gain)
    println("Estimated storage time at 100% power: $(round(storage_time, digits=2)) s")
    println()
    plot_rmt((1:length(rel_power_gain)).*4, rel_power_gain .* 100; xlabel="Time [s]", ylabel="Rel. Power Gain [%]", fig="rel_power_ref", pltctrl)
    results = JLD2.load(data_file_group_control, "results")
    print_gains(optimal_scaling)
else
    results = JLD2.load(data_file, "results")
end

if ONLINE
    println("Creating video from velocity reduction frames")
    video_path = createVideo("ff_velocity_reduction"; fps=6)
    if !isempty(video_path)
        println("✓ Created video: $video_path")
    else
        println("No velocity reduction frames found or video creation failed")
    end
end

results

