# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Main script to run a turbine group control (TGC) simulation with FLORIDyn.jl
# using a precomputed induction matrix for feed-forward control.
# TGC shall be extended to full model predictive control (MPC) in a future example

using FLORIDyn, TerminalPager, DistributedNext, DataFrames, NOMAD
if Threads.nthreads() == 1; using ControlPlots; end

settings_file = "data/2021_54T_NordseeOne.yaml"
vis_file      = "data/vis_54T.yaml"

USE_TGC = true
USE_STEP = false
USE_FEED_FORWARD = true
ONLINE = false
T_SKIP = 400 # skip first 400s of simulation for error calculation and plotting

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
sim.end_time += 420
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

function calc_induction_matrix(demand::Vector, tuning_parameters::Vector)
    # Placeholder for actual induction matrix calculation
    return zeros(size(demand, 1), size(tuning_parameters, 1))
end

function calc_axial_induction2(time, scaling::Vector; dt=DT)
    # group_id = FLORIDyn.turbine_group(ta, turbine)   
    t1 = 240.0 + dt  # Time to start increasing demand
    t2 = 960.0 + dt  # Time to reach final demand

    scaling_begin = scaling[1]
    scaling_end = scaling[2]
    # TODO calculate scaling using linear interpolation between scaling_begin and scaling_end
    scaling = scaling_begin + (scaling_end - scaling_begin) * ((time - t1) / (t2 - t1))
    scaling = max(scaling_begin, min(scaling_end, scaling))  # clamp scaling to [scaling_begin, scaling_end]
    
    demand = calc_demand(time)
    base_induction = calc_induction(demand * scaling * cp_max)

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
    
    # Apply corrections based on group
    correction = 0.0

    rel_power = calc_cp(base_induction) / cp_max + correction
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
            induction_matrix[t_idx, i + 1] = calc_axial_induction2(time, scaling)
        end
    end
    
    return induction_matrix
end

function calc_error(rel_power, demand_values, time_step)
    error = sum((rel_power[T_SKIP ÷ time_step:end] .- demand_values[T_SKIP ÷ time_step:end]).^2) / 
                length(demand_values[T_SKIP ÷ time_step:end])
    return error
end

"""
    eval_fct(x::Vector{Float64}) -> Tuple{Bool, Bool, Vector{Float64}}

Evaluation function for NOMAD optimization of the scaling parameter.

This function calculates the mean squared error between the relative power output 
and demand values for a given scaling parameter. It is designed to be minimized 
by the NOMAD optimizer.

# Arguments
- `x::Vector{Float64}`: A vector containing a single element - the scaling parameter (range: 1.0 to 2.0)

# Returns
- `success::Bool`: Always `true` to indicate successful evaluation
- `count_eval::Bool`: Always `true` to count this evaluation
- `bb_outputs::Vector{Float64}`: Vector containing the objective value (MSE)

# Global Variables Used
- `ta`: Turbine array
- `time_step`: Simulation time step
- `t_end`: End time of simulation
- `demand_values`: Target demand values
"""
function eval_fct(x::Vector{Float64})
    scaling = x  # scaling is now a vector with two elements
    
    # Calculate induction matrix with current scaling
    induction_data = calc_induction_matrix2(ta, time_step, t_end; scaling=scaling)
    
    # Run simulation and get relative power
    rel_power = run_simulation(induction_data)
    
    # Calculate error
    error = calc_error(rel_power, demand_values, time_step)
    
    bb_outputs = [error]
    success = true
    count_eval = true
    
    return (success, count_eval, bb_outputs)
end

# Set up NOMAD optimization problem
p = NomadProblem(
    2,                    # dimension (2 parameters: scaling_begin, scaling_end)
    1,                    # number of outputs (just the objective)
    ["OBJ"],             # output types: OBJ = objective to minimize
    eval_fct;            # evaluation function
    lower_bound=[1.0, 1.0],   # minimum scaling values
    upper_bound=[2.0, 2.0]    # maximum scaling values
)

# Set NOMAD options
p.options.max_bb_eval = 50      # maximum number of function evaluations
p.options.display_degree = 2    # verbosity level

# Run optimization
result = solve(p, [1.5, 1.5])  # Start from initial guess of [1.5, 1.5] 
optimal_scaling = result.x_best_feas[1:2]


induction_data = calc_induction_matrix2(ta, time_step, t_end; scaling=optimal_scaling)
rel_power = run_simulation(induction_data)
mse = calc_error(rel_power, demand_values, time_step)
println("\nRoot Mean Square Error (RMSE): $(round(sqrt(mse) * 100, digits=2))%")

plot_rmt(time_vector, [rel_power .* 100, demand_values .* 100]; xlabel="Time [s]", xlims=(T_SKIP, 1600),
         ylabel="Rel. Power Output [%]", labels=["rel_power", "rel_demand"], pltctrl)

## plot induction factor vs time for one turbine using calc_axial_induction2
# induction_factors = induction_data[:, 2]
# plot_rmt(time_vector, induction_factors; xlabel="Time [s]", ylabel="Axial Induction Factor", fig="induction", pltctrl)
