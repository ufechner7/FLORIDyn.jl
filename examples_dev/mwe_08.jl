using ControlPlots

USE_TGC = false
USE_STEP = false
USE_FEED_FORWARD = true # if false, use constant induction (no feed-forward)
ONLINE = false
T_START = 240   # time to start increasing demand
T_END   = 960   # time to reach final demand

include("../examples/calc_induction_matrix.jl")

function interpolate_scaling_lagrange(time, t1, t2, scaling::Vector{Float64})
    """Original Lagrange interpolation (can have dips)"""
    scaling_begin = scaling[1]
    scaling_mid = scaling[2]
    scaling_end = scaling[3]
    # Quadratic interpolation that satisfies:
    # time=t1 -> scaling_begin, time=(t1+t2)/2 -> scaling_mid, time=t2 -> scaling_end
    s = clamp((time - t1) / (t2 - t1), 0.0, 1.0)
    L0 = 2 * (s - 0.5) * (s - 1.0)    # basis for value at s=0
    L1 = 4 * s * (1.0 - s)            # basis for value at s=0.5
    L2 = 2 * s * (s - 0.5)            # basis for value at s=1
    result = scaling_begin * L0 + scaling_mid * L1 + scaling_end * L2
    return result
end

function interpolate_scaling(time, t1, t2, scaling::Vector{Float64}; group_id=1)
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

    """Monotonic piecewise cubic Hermite spline with C1 continuity"""
    scaling_begin = scaling[1]
    scaling_mid = scaling[2]
    scaling_end = scaling[3]
    
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
        
        result = h00 * scaling_begin + h10 * slope_begin * t_mid + 
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
        
        result = h00 * scaling_mid + h10 * slope_mid * (1.0 - t_mid) + 
                 h01 * scaling_end + h11 * slope_end * (1.0 - t_mid)
    end

    demand = calc_demand(time)
    demand_end = calc_demand(t2)
    # id_scaling interpolates between demand (1.0) and demand_end (0.0)
    # result is the time-varying scaling factor
    interpolated_demand = demand_end - (demand_end - demand) * id_scaling
    scaled_demand = result * interpolated_demand
    base_induction = calc_induction(scaled_demand * cp_max)
    
    # Debug print for a specific time point
    if abs(time - 800.0) < 2.0
        println("Time: $time, group_id: $group_id, id_scaling: $id_scaling")
        println("  demand: $demand, demand_end: $demand_end")
        println("  interpolated_demand: $interpolated_demand")
        println("  result: $result, scaled_demand: $scaled_demand")
    end

    return result, demand, scaled_demand, base_induction
end

# Example usage
time_step = 4
t_end = 1620
time_vector = 0:time_step:t_end
dt = 400
t1 = 240.0 + dt  # Time to start increasing demand
t2 = 960.0 + dt  # Time to reach final demand

# scaling = [1.1, 1.15, 1.25, 0.5, 0.8, 1.5]
scaling = [1.261, 1.285, 1.316, 0.0031, 1.994, 0]

# Calculate scaling values over time for both methods
group_id = 4
results_tuples = [interpolate_scaling(t, t1, t2, scaling; group_id=group_id) for t in time_vector]
result = [r[1] for r in results_tuples]
demand = [r[2] for r in results_tuples]
scaled_demand = [r[3] for r in results_tuples]
base_induction = [r[4] for r in results_tuples]
scaling_values_lagrange = [interpolate_scaling_lagrange(t, t1, t2, scaling) for t in time_vector]

# Plot both methods for comparison
plot_rmt(collect(time_vector), [result, demand, scaled_demand, base_induction];
         xlabel="Time [s]",
         ylabel="Scaling Factor [-]",
         title="Demand Scaling Factor vs Time",
         labels=["scaling_values_spline", "demand", "scaled_demand", "base_induction"],
         pltctrl=ControlPlots)

# Calculate base_induction for all 4 groups
base_induction_group1 = [interpolate_scaling(t, t1, t2, scaling; group_id=1)[4] for t in time_vector]
base_induction_group2 = [interpolate_scaling(t, t1, t2, scaling; group_id=2)[4] for t in time_vector]
base_induction_group3 = [interpolate_scaling(t, t1, t2, scaling; group_id=3)[4] for t in time_vector]
base_induction_group4 = [interpolate_scaling(t, t1, t2, scaling; group_id=4)[4] for t in time_vector]

# Plot base_induction vs time for all groups
plot_rmt(collect(time_vector), [base_induction_group1, base_induction_group2, base_induction_group3, base_induction_group4];
         xlabel="Time [s]",
         ylabel="Base Induction Factor [-]",
         title="Base Induction Factor vs Time by Group",
         labels=["Group 1", "Group 2", "Group 3", "Group 4"],
         fig="Base Induction by Group",
         pltctrl=ControlPlots)


