# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Main script to run a model predictive control (MPC) simulation with FLORIDyn.jl
# Currently, two modes are supported: control of all turbines with the same induction factor using CONTROL_POINTS parameters
# and group control with CONTROL_POINTS + (GROUPS-1) parameters (CONTROL_POINTS for induction corrections at different time points 
# and GROUPS-1 for individual group corrections, with the last group calculated from a constraint).

# The script uses the NOMAD.jl package for black-box optimization of the correction parameters.
# The number of groups can be 1, 2, 3, 4, 6, 8, or 12.
# The mean square error between the production and demand is minimized.

# Always run the script with GROUPS = 1 first to get a baseline result without group control. This baseline 
# is stored in the file data/mpc_puls_108_result.jld2 and used for comparison when running with group control.

# The constant MAX_STEPS can be one for debugging. With 10 or 100 you already get a rough idea of the optimization 
# progress. So far, I never needed more than about 1600 steps to converge with GROUPS = 12 and less with fewer groups.

# If you want to create a video, make sure to set ONLINE = true and run Julia with sufficient threads. On a 
# 7850X with 16 threads, the video creation needs about two hours.

# To create a bar plot, run Julia single threaded.

using Pkg
if !("NOMAD" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using FLORIDyn, TerminalPager, DistributedNext, DataFrames, NOMAD, JLD2, Statistics, Printf
using FLORIDyn: TurbineGroup, TurbineArray
if Threads.nthreads() == 1; using ControlPlots; end

const settings_file = "data/2026_108T_NordseeOne.yaml"
const vis_file      = "data/vis_108T.yaml"
const data_file               = "data/mpc_puls_108_result.jld2"
const reference_file          = "data/mpc_puls_108_reference.jld2"
const error_file              = "data/mpc_puls_108_error.jld2"
const data_file_group_control = "data/mpc_puls_108_result_group_control"

const GROUPS = 1 # for USE_HARDCODED_INITIAL_GUESS: 1, 2, 3, 4, 6, 8 or 12, otherwise any integer >= 1
CONTROL_POINTS = 7
MAX_ID_SCALING = 3.0
MAX_STEPS = 1     # maximum number black-box evaluations for NOMAD optimizer; zero means load cached results if available
USE_HARDCODED_INITIAL_GUESS = false # set to false to start from generic initial guess
USE_ADVECTION = true  
USE_PULSE = true
USE_TGC = false
USE_STEP = false
USE_FEED_FORWARD = true # if false, use constant induction (no feed-forward)
ONLINE  = false    # if true, enable online plotting during simulation and create video
const TURBULENCE = true # if true, show the added turbulence in the visualization
T_START = 240      # relative time to start increasing demand
T_END   = 2260+4000     # relative time to reach final demand
T_SHIFT = 300      # time shift the demand compared to the wind speed in seconds
T1 = 1940.0        # start time of correction spline
T2 = 5686.0        # end time of correction spline
REL_POWER = 0.90   # relative power for pulse demand
if USE_ADVECTION
    T_EXTRA = 4880+12000    # extra time in addition to sim.end_time for MPC simulation
else
    T_EXTRA = 2580+12000    # extra time in addition to sim.end_time for MPC simulation
end
MIN_INDUCTION = 0.01
MAX_DISTANCES = Float64[]
# data_file_group_control_full = data_file_group_control * '_' * string(GROUPS) * "TGs.jld2"
rel_power_ref = nothing
spline_positions = Float64[]
rel_spline_positions = Float64[]
rel_demands = Float64[]
max_powers = Float64[]
scaled_demands = Float64[]

const GROUP_CONTROL = (GROUPS != 1)
if USE_HARDCODED_INITIAL_GUESS
    @assert(GROUPS in (1, 2, 3, 4, 6, 8, 12), "GROUPS must be 1, 2, 3, 4, 6, 8, or 12")
else
    @assert(GROUPS >= 1, "GROUPS must be at least 1")
end
if MAX_STEPS == 0
   const SIMULATE = false # if false, load cached results if available
else
   const SIMULATE = true
end
if TURBULENCE
    msr = AddedTurbulence
else
    msr = VelReduction
end
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
md::DataFrame = DataFrame()

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
induction = 1/3  # Use Betz limit (optimal induction factor)
set_induction!(ta, induction)

time_step = sim.time_step  # seconds
t_end = sim.end_time - sim.start_time  # relative end time in seconds

function calc_vel(vis, ta::TurbineArray, start_time, t_end)
    n_turbines = size(ta.pos, 1)
    time_vector = start_time:time_step:t_end
    u0 = calc_wind.(Ref(vis), time_vector.-start_time)
    c_true = mean(u0)     # m/s
    # `wind_data::Matrix{Float64}`: Matrix where each row is time, U_T0, U_T1, ... U_Tn.
    wind_data = zeros(Float64, length(time_vector), n_turbines + 1)
    wind_data[:, 1] = collect(time_vector)
    for i in eachindex(ta.pos[:, 1])
        # calculate the x position of each turbine
        x_pos = ta.pos[i, 1]
        delay_steps = round(Int, x_pos / c_true / time_step)  # Number of time steps for delay
        # Create delayed signal: u_x[i] = u0[i - delay_steps]
        u_x = zeros(Float64, length(u0))
        for j in eachindex(u0)
            src_idx = max(1, j - delay_steps)
            u_x[j] = u0[src_idx]
        end
        wind_data[:, i + 1] .= u_x
    end
    return wind_data
end

"""
    u_mean(wind_data) -> Vector{Float64}

Calculate the cubic mean wind speed over all turbines for each time step.

# Arguments
- `wind_data::Matrix{Float64}`: Matrix where column 1 is time and columns 2:end are wind speeds for each turbine

# Returns
- `Vector{Float64}`: Cubic mean wind speed for each time step: (mean(u^3))^(1/3)
"""
function u_mean(wind_data)
    # Extract wind speeds (skip first column which is time)
    wind_speeds = wind_data[:, 2:end]
    n_turbines = size(wind_speeds, 2)
    
    # Calculate cubic mean for each time step (row)
    return [cbrt(sum(wind_speeds[i, :] .^ 3) / n_turbines) for i in axes(wind_speeds, 1)]
end

"""
    u_mean(wind_data, time) -> Float64

Calculate the cubic mean wind speed over all turbines at a specific time point.

# Arguments
- `wind_data::Matrix{Float64}`: Matrix where column 1 is time and columns 2:end are wind speeds for each turbine
- `time::Float64`: Time point at which to calculate the mean wind speed

# Returns
- `Float64`: Cubic mean wind speed at the specified time: (mean(u^3))^(1/3)
"""
function u_mean(wind_data, time)
    # Find the row index corresponding to the given time
    time_vec = wind_data[:, 1]
    idx = findfirst(t -> t >= time, time_vec)
    
    if isnothing(idx)
        idx = length(time_vec)  # Use last time point if time is beyond range
    end
    
    # Extract wind speeds for all turbines at this time point
    wind_speeds = wind_data[idx, 2:end]
    n_turbines = length(wind_speeds)
    
    # Calculate cubic mean: (mean(u^3))^(1/3)
    return cbrt(sum(wind_speeds .^ 3) / n_turbines)
end


# Calculate demand for each time point
time_vector = 0:time_step:t_end
wind_data = [calc_wind(vis, t) for t in time_vector]
if USE_ADVECTION
    # Set up wind velocity interpolation BEFORE creating induction matrix and settings
    wind.input_vel = "InterpTurbine"
    wind.vel = calc_vel(vis, ta, sim.start_time, sim.end_time)
else
    # Set up wind velocity interpolation BEFORE creating induction matrix and settings
    wind.input_vel = "Interpolation"
    wind.vel = calc_vel(vis, sim.start_time, sim.end_time)
end
global demand_ref, demand_abs
if USE_ADVECTION && !isfile(reference_file)
    error("USE_ADVECTION is true but reference file '$reference_file' not found. " *
          "Set USE_ADVECTION = false and run first to generate reference data.")
end

if isfile(reference_file) && USE_ADVECTION
    # Load reference relative power for error calculation
    ref_data = JLD2.jldopen(reference_file, "r") do file
        Dict(
            "time_vector" => file["results"]["time_vector"],
            "total_power" => file["results"]["total_power"]
        )
    end
    demand_ref = ref_data["total_power"]
    # Replace values below 100 MW with last valid value
    let last_valid = 0.0, n_replaced = 0
        for i in 1:length(demand_ref)
            if demand_ref[i] >= 100 && ! isnan(demand_ref[i])
                last_valid = demand_ref[i]
            elseif last_valid > 0
                demand_ref[i] = last_valid
                n_replaced += 1
            end
        end
        if n_replaced > 0
            println("Replaced $n_replaced values below 100 MW with last valid value")
        end
    end
       
    println("Loaded demand_ref with $(length(demand_ref)) time steps")
    println("demand_ref range: $(minimum(demand_ref)) to $(maximum(demand_ref)) MW")
    
    # apply T_SHIFT to demand_abs
    if T_SHIFT != 0
        n_shift_steps = round(Int, T_SHIFT / time_step)
        demand_abs = vcat(zeros(Float64, n_shift_steps), demand_ref[1:end - n_shift_steps])
    else
        demand_abs = copy(demand_ref)
    end
    demand_abs .*= REL_POWER
else
    # Initialize demand_ref as empty for non-advection cases
    demand_ref = Float64[]
    # Calculate demand in absolute power (Watts)
    demand_abs = [calc_demand(vis, t; t_shift=T_SHIFT, rel_power=REL_POWER) for t in time_vector]
end

# Replace values below 100 MW with last valid value
let last_valid = 0.0, n_replaced = 0
    for i in 1:length(demand_abs)
        if demand_abs[i] >= 100
            last_valid = demand_abs[i]
        elseif last_valid > 0
            demand_abs[i] = last_valid
            n_replaced += 1
        end
    end
    if n_replaced > 0
        println("Replaced $n_replaced values below 100 MW with last valid value")
    end
end

demand_data = demand_abs 

# For initial setup, use calc_induction_matrix (only affects pre-optimization visualization)
# During optimization, calc_induction_matrix2 will be used with proper group handling
# Do NOT call prepareSimulation here - it will be called in run_simulation()
# Create constant induction matrix with value 1/3 (Betz limit)
time_vec = 0:time_step:t_end
n_time_steps = length(time_vec)
n_turbines = size(ta.pos, 1)
con.induction_data = zeros(Float64, n_time_steps, n_turbines + 1)
con.induction_data[:, 1] = collect(time_vec)
con.induction_data[:, 2:end] .= 1/3  # Set constant induction of 1/3 for all turbines

# create settings struct with automatic parallel/threading detection
set = Settings(wind, sim, con, Threads.nthreads() > 1, Threads.nthreads() > 1)
if USE_FEED_FORWARD
    set.induction_mode = Induction_TGC()
else
    set.induction_mode = Induction_Constant()
end

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
    return max_power
end

# This function implements the "model" in the block diagram.
function run_simulation(set_induction::AbstractMatrix; enable_online=false, msr=msr)
    global set, wind, con, floridyn, floris, sim, ta, vis
    con.induction_data = set_induction
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    # Only enable online visualization if explicitly requested (to avoid NaN issues during optimization)
    vis.online = enable_online
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris; msr=msr)
    
    # Calculate total wind farm power by grouping by time and summing turbine powers
    total_power_df = combine(groupby(md, :Time), :PowerGen => sum => :TotalPower)
    
    # Return absolute power in Watts
    abs_power = total_power_df.TotalPower
    
    return abs_power
end

# Separate function to run simulation and return full DataFrame for plotting
function run_simulation_full(set_induction::AbstractMatrix; enable_online=false, msr=msr)
    global set, wind, con, floridyn, floris, sim, ta, vis
    con.induction_data = set_induction
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    vis.online = enable_online
    wf, md, mi = run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris; msr=msr)
    return md
end


"""
    calc_axial_induction2(vis, time, correction::Vector; group_id=nothing) -> (corrected_induction, distance)

Calculate the axial induction factor for a turbine using optimizable correction parameters.

This function computes the axial induction factor based on:
1. Time-dependent correction via cubic Hermite spline interpolation of control points
2. Optional group-specific correction factors for individual turbine group control
3. Demand-based adjustment with correction for power coefficient nonlinearity

# Arguments
- `vis`: [`Vis`](@ref) object containing visualization settings (uses `t_skip`)
- `time::Float64`: Current simulation time in seconds
- `correction::Vector{Float64}`: Optimization parameters vector containing:
  - Elements 1 to CONTROL_POINTS: Time-dependent correction control points
  - Elements CONTROL_POINTS+1 to end: Group-specific correction factors (if GROUP_CONTROL)
- `group_id::Union{Int,Nothing}`: Turbine group identifier (1 to GROUPS), or `nothing` for no group control

# Returns
- `corrected_induction::Float64`: Computed axial induction factor, clamped to [MIN_INDUCTION, BETZ_INDUCTION]
- `distance::Float64`: Constraint violation distance (positive if scaled_demand > 1.0, else 0.0)

# Details
The function operates in several stages:
1. Extracts group-specific correction factor `id_correction` from the correction vector (if applicable)
   - For groups 1 to GROUPS-1: directly from `correction[CONTROL_POINTS + group_id]`
   - For the last group (GROUPS): calculated as `GROUPS * MAX_ID_SCALING / 2.0 - sum(correction[(CONTROL_POINTS+1):end])`
     to reduce the number of optimization variables by one
2. Computes normalized time parameter `s` ∈ [0,1] between T_START and T_END
3. Interpolates time-dependent correction using [`interpolate_hermite_spline`](@ref)
4. Adjusts demand by group-specific correction and applies time-dependent correction
5. Converts scaled demand to induction, applies power coefficient correction
6. Ensures minimum induction to avoid numerical issues in wake model

# Global Constants Used
- `CONTROL_POINTS`: Number of time-dependent control points
- `GROUPS`: Number of turbine groups
- `MAX_ID_SCALING`: Maximum allowed group correction factor
- `T_START`: Time offset to start ramping demand (relative to `vis.t_skip`)
- `T_END`: Time offset to reach final demand (relative to `vis.t_skip`)
- `MIN_INDUCTION`: Minimum induction to prevent NaN in FLORIS wake model
- `BETZ_INDUCTION`: Maximum theoretical induction (Betz limit)

# See Also
- [`interpolate_hermite_spline`](@ref): Performs cubic Hermite spline interpolation
- [`calc_induction_matrix2`](@ref): Uses this function to build induction matrices
"""
function calc_axial_induction2(vis, time, correction::Vector; group_id=nothing)
    global spline_positions, rel_spline_positions, rel_demands
    distance = 0.0
    id_correction = 1.0
    if length(correction) > CONTROL_POINTS && !isnothing(group_id) && group_id >= 1
        if group_id <= GROUPS - 1
            id_correction = correction[CONTROL_POINTS + group_id]
        elseif group_id == GROUPS
            # Last group: calculate as GROUPS * MAX_ID_SCALING / 2.0 minus sum of groups 1 to GROUPS-1
            id_correction = GROUPS * MAX_ID_SCALING / 2.0 - sum(correction[(CONTROL_POINTS+1):end])
        end
        id_correction = clamp(id_correction, 0.0, MAX_ID_SCALING)
    end

    t1 = T1
    t2 = T2
    
    # Calculate s_positions equally spaced between 0 and 1
    s_positions = range(0.0, 1.0, length=CONTROL_POINTS) |> collect
    rel_spline_positions = s_positions
    
    # Calculate absolute time positions for spline control points
    spline_positions = [t1 + s * (t2 - t1) for s in s_positions]
    
    # Calculate normalized time parameter s for interpolation
    # Clamp time for s calculation, but preserve original time for demand/wind
    time_clamped = max(time, t1)
    s = clamp((time_clamped - t1) / (t2 - t1), 0.0, 1.0)
    # Perform piecewise cubic Hermite spline interpolation
    correction_result = interpolate_hermite_spline(s, correction[1:CONTROL_POINTS], s_positions)
    # correction_result = 1.0
    
    demand = calc_demand(vis, time; t_shift=T_SHIFT, rel_power=REL_POWER)
    scaled_demand = correction_result * demand
    max_power = calc_demand(vis, time; t_shift=0.0, rel_power=1.0)
    
    # Protect against division by very small max_power values
    if max_power < 1e6  # Less than 1 MW
        rel_demand = 0.0
    else
        rel_demand = scaled_demand / max_power
    end
    
    # Apply group-specific correction
    rel_demand *= id_correction
    if rel_demand > 1.0
        rel_demand = 1.0
    end
    push!(rel_demands, rel_demand)
    push!(max_powers, max_power)
    push!(scaled_demands, scaled_demand)
    base_induction = calc_induction(rel_demand * cp_max)

    rel_power = calc_cp(base_induction) / cp_max
    corrected_induction = calc_induction(rel_power * cp_max)
    
    # Ensure minimum induction to avoid numerical issues in FLORIS (NaN from zero induction)
    # Minimum value of MIN_INDUCTION ensures the wake model has valid inputs
    corrected_induction = max(MIN_INDUCTION, min(BETZ_INDUCTION, corrected_induction))
    # if corrected_induction >= 0.28
    #     @warn "Corrected induction is very high: $corrected_induction, $rel_demand, $scaled_demand, $max_power at time $time s"
    # end
    return corrected_induction, distance
end

function calc_induction_matrix2(vis, ta, time_step, t_end; correction)
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
        for i in eachindex(ta.pos[:, 1])
            group_id = FLORIDyn.turbine_group(ta, i)
            time_shifted = time - time_step  # shift back by one time step to suppress spike
            axial_induction, distance = calc_axial_induction2(vis, time_shifted, correction; group_id=group_id)
            induction_matrix[t_idx, i+1] = axial_induction
            max_distance = max(max_distance, distance)
        end
    end

    return induction_matrix, max_distance
end

function eval_fct(x::Vector{Float64})
    global rel_power
    correction = x  # correction is now a vector with parameters
    print(".")  # progress indicator
    
    # Calculate induction matrix with current correction
    induction_data, max_distance = calc_induction_matrix2(vis, ta, time_step, t_end; correction=correction)
    if max_distance > 0.0
        push!(MAX_DISTANCES, max_distance)
    end

    # Run simulation and get absolute power
    abs_power = run_simulation(induction_data)
    
    # Calculate error
    error = calc_error(vis, abs_power, demand_data, time_step)
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
    n_total_params = CONTROL_POINTS + n_group_params  # CONTROL_POINTS global correction + (GROUPS-1) group correction
    
    # Create lower and upper bounds dynamically
    lower_bound = vcat(fill(0.8, CONTROL_POINTS), fill(0.0, n_group_params))
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
        CONTROL_POINTS,      # dimension (CONTROL_POINTS parameters: correction at CONTROL_POINTS time points)
        1,                   # number of outputs (just the objective)
        ["OBJ"],             # output types: OBJ = objective to minimize
        eval_fct;            # evaluation function
        lower_bound=fill(0.8, CONTROL_POINTS),   # minimum correction values
        upper_bound=vcat([1.2], fill(1.2, CONTROL_POINTS - 1))    # maximum correction values
    )

    # Set NOMAD options
    p.options.max_bb_eval = MAX_STEPS      # maximum number of function evaluations
    p.options.display_degree = 2           # verbosity level
end


function calc_error(vis, abs_power, demand_data, time_step)
    # Start index after skipping initial transient; +1 because Julia is 1-based
    i0 = Int(floor(vis.t_skip / time_step)) + 1
    # Clamp to valid range
    i0 = max(1, i0)
    n = min(length(abs_power), length(demand_data)) - i0 + 1
    if n <= 0
        error("calc_error: empty overlap after skip; check vis.t_skip and lengths (abs_power=$(length(abs_power)), demand=$(length(demand_data)), i0=$(i0))")
    end
    p = @view abs_power[i0:i0 + n - 1]
    d = @view demand_data[i0:i0 + n - 1]
    return sqrt(sum((p .- d) .^ 2) / length(d))
end

include("mpc_plotting.jl")

# Prepare simulation to get wf and floris (needed for calc_axial_induction2)
wf, wind_prep, sim_prep, con_prep, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)

# demand_data is already in absolute power (Watts), no conversion needed

if SIMULATE
    println("Starting NOMAD optimization with max $(p.options.max_bb_eval) evaluations...")
    if GROUP_CONTROL
        if GROUPS == 2
            # Hardcoded initial guess from previous runs (originally for 11 control points)
            x0_full =  [1.0, 1.1, 1.1, 1.1, 1.3, 1.1, 1.6, 2.0, 1.5, 1.502, 1.4, 1.05]
            # Adjust to current CONTROL_POINTS: take first CONTROL_POINTS + (GROUPS-1) elements
            x0 = x0_full[1:(CONTROL_POINTS + GROUPS - 1)]
        elseif GROUPS == 6
            # Hardcoded initial guess from previous runs
            x0 = [0.94518303918679, 0.93269599433136, 1.00092898512453, 1.01904687495053, 1.47398102169963, 1.99998018368223, 1.07802139832096, 1.2224803349927, 0.82466261195206, 0.98282706090237, 0.8370533275206, 0.85516767133612]
        else
            # For group control, use generic initial guess with CONTROL_POINTS corrections + (GROUPS-1) group scalings
            x0 = vcat(fill(1.0, CONTROL_POINTS), fill(1.0, GROUPS - 1))
        end
    else
        x0 = [0.99055661760185, 0.95955086565182, 0.98819178660071, 0.9856740943524, 1.03187884206546, 1.00458036586437, 0.99286704613765]
    end
    result = solve(p, x0)
    optimal_correction = result.x_best_feas
    println("\nNOMAD optimization completed.")
    println("Best correction: ", optimal_correction)
    induction_data, max_distance = calc_induction_matrix2(vis, ta, time_step, t_end; correction=optimal_correction)
    con.induction_data = induction_data  # Update controller with optimized induction data
else
    # Load cached results when not simulating
    println("Loading cached optimization results...")
    if GROUP_CONTROL
        data_file_group_control_full = data_file_group_control * "_" * string(GROUPS) * "TGs.jld2"
        if isfile(data_file_group_control_full)
            optimal_correction = JLD2.load(data_file_group_control_full, "optimal_correction")
            println("Loaded optimal correction from $data_file_group_control_full")
        else
            @warn "Cache file not found: $data_file_group_control_full. Using generic correction."
            n_group_params = GROUPS - 1
            optimal_correction = vcat(fill(1.0, CONTROL_POINTS), fill(1.0, n_group_params))
        end
    else
        if isfile(data_file)
            optimal_correction = JLD2.load(data_file, "optimal_correction")
            println("Loaded optimal correction from $data_file")
        else
            @warn "Cache file not found: $data_file. Using generic correction."
            optimal_correction = ones(CONTROL_POINTS)
        end
    end
    induction_data, max_distance = calc_induction_matrix2(vis, ta, time_step, t_end; correction=optimal_correction)
    con.induction_data = induction_data  # Update controller with optimized induction data
end

# Verify induction_data before running simulation
println("Induction data size: $(size(con.induction_data))")
if any(isnan.(con.induction_data))
    error("con.induction_data contains NaN values! Check calc_induction_matrix2 or wind data.")
end

# Check wind velocity data for NaN values
if any(isnan.(wind.vel))
    @warn "wind.vel contains $(sum(isnan.(wind.vel))) NaN values!"
    nan_rows = findall(any(isnan.(wind.vel), dims=2)[:])
    println("NaN found in wind.vel at rows: $nan_rows")
    println("Time values at NaN rows: $(wind.vel[nan_rows, 1])")
end

# Check time points around problematic area before simulation
problematic_time = 28236.0
time_idx = findfirst(t -> abs(t - problematic_time) < time_step, con.induction_data[:, 1])
if !isnothing(time_idx)
    println("\nDiagnosing time around $problematic_time seconds (index $time_idx):")
    start_idx = max(1, time_idx - 2)
    end_idx = min(size(con.induction_data, 1), time_idx + 2)
    println("Induction data rows $start_idx to $end_idx:")
    println("Time: ", con.induction_data[start_idx:end_idx, 1])
    println("Sample induction values (turbine 1): ", con.induction_data[start_idx:end_idx, 2])
    
    # Check wind velocity at this time
    if USE_ADVECTION
        wind_idx = findfirst(t -> abs(t - problematic_time) < time_step, wind.vel[:, 1])
        if !isnothing(wind_idx)
            wind_start = max(1, wind_idx - 2)
            wind_end = min(size(wind.vel, 1), wind_idx + 2)
            println("Wind vel rows $wind_start to $wind_end:")
            println("Time: ", wind.vel[wind_start:wind_end, 1])
            println("Sample wind speeds (turbine 1): ", wind.vel[wind_start:wind_end, 2])
        end
    end
end

# Run final simulation and get full DataFrame for plotting
md = run_simulation_full(con.induction_data)

# Check if md.Time contains NaN and diagnose
if any(isnan.(md.Time))
    nan_indices = findall(isnan.(md.Time))
    println("\n⚠️  WARNING: md.Time contains $(length(nan_indices)) NaN value(s)!")
    println("NaN found at row indices: $(nan_indices[1:min(10, length(nan_indices))])")
    
    # Show surrounding rows for context
    if length(nan_indices) > 0
        idx = nan_indices[1]
        start_idx = max(1, idx - 2)
        end_idx = min(nrow(md), idx + 2)
        println("\nContext around first NaN (rows $start_idx to $end_idx):")
        println(md[start_idx:end_idx, [:Time, :PowerGen]])
        
        # Try to identify which turbines have issues
        turbine_col = :Turbine in names(md) ? :Turbine : :turbineID
        if turbine_col in names(md)
            nan_turbines = unique(md[nan_indices, turbine_col])
            println("\nTurbines with NaN: $(nan_turbines[1:min(10, length(nan_turbines))])")
        end
    end
    
    # Filter out NaN rows for plotting
    println("\nFiltering out $(length(nan_indices)) NaN rows from md.Time for plotting...")
    md = md[.!isnan.(md.Time), :]
    println("Remaining rows: $(nrow(md))")
end

# Plot wind speed vs time (using relative time)
if USE_ADVECTION
    # Extract time from wind.vel and compute relative time
    wind_time_rel = wind.vel[:, 1] .- sim.start_time
    u_mean_vals = u_mean(wind.vel)
    
    # Ensure all arrays have the same length by finding minimum
    min_len = min(length(time_vector), length(wind_data), length(u_mean_vals))
    println("Plotting lengths: time_vector=$(length(time_vector)), wind_data=$(length(wind_data)), u_mean=$(length(u_mean_vals)), using min_len=$min_len")
    
    plot_rmt(collect(time_vector)[1:min_len], [wind_data[1:min_len], u_mean_vals[1:min_len]]; 
        xlabel="Time [s]", xlims=(vis.t_skip, collect(time_vector)[min_len]),
        ylabel="v_wind [m/s]", labels=["u_inflow","u_mean"], fig="v_wind", title="Wind speed vs time", pltctrl)
else
    plot_rmt(collect(time_vector), wind_data; xlabel="Time [s]", xlims=(vis.t_skip, time_vector[end]),
        ylabel="v_wind [m/s]", fig="v_wind", title="Wind speed vs time", pltctrl)
end

# Plot total power output and demand vs time
if "PowerGen" in names(md)
    # Group by time and sum power across all turbines
    time_points = unique(md.Time)
    total_power = [sum(md[md.Time .== t, "PowerGen"]) for t in time_points]
    
    # Convert absolute time to relative time for consistent x-axis
    time_points_rel = time_points .- sim.start_time

    if GROUPS==1 && MAX_STEPS==0
        # store reference power
        JLD2.jldsave(reference_file; results=Dict(
            "time_vector" => collect(time_points_rel),
            "total_power" => total_power)
        )
        @info "Reference results saved to $reference_file"
    end
    
    # Calculate demand for all time points - ensure matching lengths
    min_power_len = min(length(time_points_rel), length(total_power), length(demand_abs))
    demand_power = demand_abs[1:min_power_len]
    
    plot_rmt(collect(time_points_rel)[1:min_power_len], [total_power[1:min_power_len], demand_power]; 
        xlabel="Time [s]", xlims=(vis.t_skip, time_points_rel[min_power_len]),
        ylabel="Total Power [MW]", fig="total_power", title="Total power output and demand vs time", 
        labels=["Power Output", "Demand"], pltctrl)
end

# Plot axial induction vs time using calc_axial_induction2
# Create a simple correction vector with no correction (all 1.0)
if SIMULATE
    correction = optimal_correction
else
    if GROUP_CONTROL
        n_group_params = GROUPS - 1
        correction = vcat(fill(1.0, CONTROL_POINTS), fill(1.0, n_group_params))
    else
        correction = ones(CONTROL_POINTS)
    end
end

if GROUP_CONTROL && GROUPS > 1
    local induction_values = Float64[]
    local induction
    # Plot induction for each turbine group
    group_inductions = [Float64[] for _ in 1:GROUPS]
    group_labels = ["Group $group_id" for group_id in 1:GROUPS]
    
    for group_id in 1:GROUPS
        for t in time_vector
            induction, _ = calc_axial_induction2(vis, t, correction; group_id=group_id)
            push!(group_inductions[group_id], induction)
        end
    end
    
    plot_rmt(collect(time_vector), group_inductions;
        xlabel="Time [s]",
        ylabel="Axial Induction Factor [-]",
        title="Axial induction factor vs time",
        labels=group_labels,
        fig="axial_induction",
        pltctrl=pltctrl)
else
    # Plot induction for all turbines (using group_id=1 or nothing)
    induction_values = Float64[]
    for t in time_vector
        induction, _ = calc_axial_induction2(vis, t, correction; group_id=1)
        push!(induction_values, induction)
    end
    plot_rmt(collect(time_vector), induction_values; xlabel="Time [s]", xlims=(vis.t_skip, time_vector[end]),
        ylabel="Axial Induction Factor [-]", fig="axial_induction", title="Axial induction factor vs time", pltctrl)
end

plot_correction_curve(correction, rel_spline_positions; t1=T1, t2=T2)

# Calculate absolute power from the final simulation for error calculation
total_power_df = combine(groupby(md, :Time), :PowerGen => sum => :TotalPower)
abs_power = total_power_df.TotalPower

error = calc_error(vis, abs_power, demand_data, time_step)
mean_demand = mean(demand_data[Int(floor(vis.t_skip / time_step)) + 1:end])
println("\nFinal error after optimization: $(round(error/mean_demand * 100, digits=2)) %")
println("\nFinal error after optimization: $(round(error, sigdigits=4)) MW")

# # Test case for calc_induction_matrix2: plot induction matrix per turbine group
# println("\nTesting calc_induction_matrix2...")
# test_correction = ones(CONTROL_POINTS)  # Use unity correction for testing
# induction_matrix, max_distance = calc_induction_matrix2(vis, ta, time_step, t_end; correction=test_correction)
# println("Generated induction matrix: $(size(induction_matrix)) (time_steps × turbines+1)")
# println("Max constraint distance: $max_distance")

# # Extract time from first column
# matrix_time = induction_matrix[:, 1]

# if GROUP_CONTROL && GROUPS > 1
#     # Plot average induction per group from the matrix
#     group_matrix_inductions = [Float64[] for _ in 1:GROUPS]
#     group_matrix_labels = ["Group $group_id" for group_id in 1:GROUPS]
    
#     for group_id in 1:GROUPS
#         # Find turbines in this group
#         turbine_indices = [i for i in 1:n_turbines if FLORIDyn.turbine_group(ta, i) == group_id]
#         if !isempty(turbine_indices)
#             # Average induction across turbines in this group (columns are turbine_index + 1)
#             avg_induction = mean(induction_matrix[:, turbine_indices .+ 1], dims=2)[:]
#             group_matrix_inductions[group_id] = avg_induction
#         else
#             group_matrix_inductions[group_id] = zeros(size(induction_matrix, 1))
#         end
#     end
    
#     plot_rmt(collect(matrix_time), group_matrix_inductions;
#         xlabel="Time [s]",
#         ylabel="Axial Induction Factor [-]",
#         title="Induction matrix test (calc_induction_matrix2)",
#         labels=group_matrix_labels,
#         fig="induction_matrix_test",
#         pltctrl=pltctrl)
# else
#     # Plot average induction across all turbines
#     avg_matrix_induction = mean(induction_matrix[:, 2:end], dims=2)[:]
#     plot_rmt(collect(matrix_time), avg_matrix_induction; xlabel="Time [s]", xlims=(vis.t_skip, matrix_time[end]),
#         ylabel="Axial Induction Factor [-]", fig="induction_matrix_test", 
#         title="Induction matrix test (calc_induction_matrix2)", pltctrl)
# end
