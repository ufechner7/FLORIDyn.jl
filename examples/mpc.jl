# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Main script to run a model predictive control (MPC) simulation with FLORIDyn.jl
# Currently, two modes are supported: control of all turbines with the same induction factor using 3 parameter
# and group control with 6 parameters (3 for induction scaling and 3 for individual group scaling).
# The mean square error between the production and demand is minimized.

using Pkg
if ! ("NOMAD" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using FLORIDyn, TerminalPager, DistributedNext, DataFrames, NOMAD, JLD2
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"
data_file               = "data/mpc_result.jld2"
data_file_group_control = "data/mpc_result_group_control.jld2"

GROUP_CONTROL = true  # if false, use 3-parameter control for all turbines; if true, use 6-parameter group control
SIMULATE = false      # if false, load cached results if available
MAX_STEPS = 100       # maximum number black-box evaluations for NOMAD optimizer
USE_TGC = false
USE_STEP = false
USE_FEED_FORWARD = true # if false, use constant induction (no feed-forward)
ONLINE = false
T_SKIP = 400    # skip first 400s of simulation for error calculation and plotting
T_START = 240   # time to start increasing demand
T_END   = 960   # time to reach final demand
T_EXTRA = 1520  # extra time in addition to sim.end_time for MPC simulation

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
    id_scaling = 1.0
    if length(scaling) > 3
        if group_id == 1
            id_scaling = scaling[4]
        elseif group_id == 2
            id_scaling = scaling[5]
        elseif group_id == 3
            id_scaling = scaling[6]
        else
            id_scaling = 4.0 - (scaling[4] + scaling[5] + scaling[6])
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
    
    # Monotonic piecewise cubic Hermite spline interpolation with C1 continuity
    # Normalized time parameter
    s = clamp((time - t1) / (t2 - t1), 0.0, 1.0)
    
    # Split into two segments at s=0.5
    t_mid = 0.5
    
    # Calculate slopes at each point using finite differences
    slope1 = 2 * (scaling_mid - scaling_begin)  # slope from begin to mid
    slope2 = 2 * (scaling_end - scaling_mid)    # slope from mid to end
    
    # Derivative at beginning (use slope of first segment)
    slope_begin = slope1
    
    # Derivative at midpoint (average of adjacent slopes for C1 continuity)
    slope_mid = (slope1 + slope2) / 2.0
    
    # Derivative at end (use slope of second segment)
    slope_end = slope2
    
    if s <= t_mid
        # First segment: [0, 0.5]
        s_local = s / t_mid  # normalize to [0, 1]
        # Hermite interpolation: f(0)=scaling_begin, f(1)=scaling_mid
        # f'(0)=slope_begin, f'(1)=slope_mid
        h00 = 2*s_local^3 - 3*s_local^2 + 1
        h10 = s_local^3 - 2*s_local^2 + s_local
        h01 = -2*s_local^3 + 3*s_local^2
        h11 = s_local^3 - s_local^2
        
        scaling_result = h00 * scaling_begin + h10 * slope_begin * t_mid + 
                        h01 * scaling_mid + h11 * slope_mid * t_mid
    else
        # Second segment: [0.5, 1.0]
        s_local = (s - t_mid) / (1.0 - t_mid)  # normalize to [0, 1]
        # Hermite interpolation: f(0)=scaling_mid, f(1)=scaling_end
        # f'(0)=slope_mid, f'(1)=slope_end
        h00 = 2*s_local^3 - 3*s_local^2 + 1
        h10 = s_local^3 - 2*s_local^2 + s_local
        h01 = -2*s_local^3 + 3*s_local^2
        h11 = s_local^3 - s_local^2
        
        scaling_result = h00 * scaling_mid + h10 * slope_mid * (1.0 - t_mid) + 
                        h01 * scaling_end + h11 * slope_end * (1.0 - t_mid)
    end
    
    demand = calc_demand(time)
    demand_end = calc_demand(t2)
    interpolated_demand = demand_end - (demand_end - demand) * id_scaling
    scaled_demand = scaling_result * interpolated_demand
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
    return max(0.0, min(BETZ_INDUCTION, corrected_induction))
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
    for (t_idx, time) in enumerate(time_vector)
        for i in 1:n_turbines
            group_id = FLORIDyn.turbine_group(ta, i)
            induction_matrix[t_idx, i + 1] = calc_axial_induction2(time, scaling; group_id=group_id)
        end
    end
    
    return induction_matrix
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
    induction_data = calc_induction_matrix2(ta, time_step, t_end; scaling=scaling)
    
    # Run simulation and get relative power
    rel_power = run_simulation(induction_data)
    
    # Calculate error
    error = calc_error(rel_power, demand_values, time_step)
    
    # Add constraint if GROUP_CONTROL is true
    if GROUP_CONTROL
        # Constraint: x[4] + x[5] + x[6] <= 4
        # For NOMAD, constraints should be <= 0, so we formulate as:
        # x[4] + x[5] + x[6] - 4 <= 0
        constraint = x[4] + x[5] + x[6] - 4.0
        bb_outputs = [error, constraint]
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
        6,                    # dimension (6 parameters: scaling_begin, scaling_mid, scaling_end, id_scaling)
        2,                    # number of outputs (objective + 1 constraint)
        ["OBJ", "PB"],       # output types: OBJ = objective to minimize, PB = progressive barrier constraint
        eval_fct;            # evaluation function
        lower_bound=[1.0, 1.0, 1.0, 0.0, 0.0, 0.0],   # minimum scaling values
        upper_bound=[2.0, 2.0, 2.0, 2.0, 2.0, 2.0]    # maximum scaling values
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
        result = solve(p, [1.261, 1.285, 1.316, 0.0031, 1.994, 0])
        results_ref = JLD2.load(data_file, "results") 
        optimal_scaling = result.x_best_feas[1:6]
    else
        result = solve(p, [1.5, 1.5, 1.5])  # Start from initial guess of [1.5, 1.5, 1.5]
        optimal_scaling = result.x_best_feas[1:3]
    end

    induction_data = calc_induction_matrix2(ta, time_step, t_end; scaling=optimal_scaling)
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

    # Prepare containers for four groups
    group_data = [Float64[] for _ in 1:4]

    # Since all turbines in a group have identical induction, pick one representative turbine per group
    group_indices = [findfirst(i -> FLORIDyn.turbine_group(ta, i) == g, 1:n_turbines) for g in 1:4]
    for g in 1:4
        idx = group_indices[g]
        if isnothing(idx)
            group_data[g] = fill(0.0, n_time_steps)
        else
            group_data[g] = induction_data[:, idx + 1]
        end
    end

    group_labels = ["Group 1", "Group 2", "Group 3", "Group 4"]
    plot_rmt(time_vec_ind, group_data;
             xlabel="Time [s]",
             ylabel="Axial Induction Factor [-]",
             title="Average Axial Induction Factor vs Time by Turbine Group",
             labels=group_labels,
             fig="Induction by Group",
             pltctrl=pltctrl)
end

if GROUP_CONTROL
    # calculate rel_power-rel_power_ref
    start_index = Int(floor((T_SKIP+T_START+(T_END-T_SKIP-T_START)/2) / time_step)) + 1
    rel_power_gain = rel_power[start_index:end-1] .- rel_power_ref[start_index:end]
    plot_rmt((1:length(rel_power_gain)).*4, rel_power_gain .* 100; xlabel="Time [s]", ylabel="Rel. Power Gain [%]", fig="rel_power_ref", pltctrl)
    results = JLD2.load(data_file_group_control, "results")
else
    results = JLD2.load(data_file, "results")
end

