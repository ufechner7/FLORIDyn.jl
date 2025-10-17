using ControlPlots

function interpolate_scaling_lagrange(time::Float64, t1::Float64, t2::Float64, scaling::Vector{Float64})
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

function interpolate_scaling(time::Float64, t1::Float64, t2::Float64, scaling::Vector{Float64})
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
    
    return result
end

# Example usage
t1 = 0.0
t2 = 10.0
scaling = [1.0, 1.2, 2.0]

# Calculate scaling values over time for both methods
time_vector = 0.0:0.1:10.0
scaling_values_spline = [interpolate_scaling(t, t1, t2, scaling) for t in time_vector]
scaling_values_lagrange = [interpolate_scaling_lagrange(t, t1, t2, scaling) for t in time_vector]

# Print some values
println("\nComparison of interpolation methods:")
println("Time\tSpline\t\tLagrange")
for time in 0.0:1.0:10.0
    s_spline = interpolate_scaling(time, t1, t2, scaling)
    s_lagrange = interpolate_scaling_lagrange(time, t1, t2, scaling)
    println("$time\t$(round(s_spline, digits=4))\t\t$(round(s_lagrange, digits=4))")
end

# Plot both methods for comparison
plot_rmt(collect(time_vector), [scaling_values_spline, scaling_values_lagrange];
         xlabel="Time [s]",
         ylabel="Scaling Factor [-]",
         title="Demand Scaling Factor vs Time",
         labels=["Quadratic Spline", "Lagrange"],
         pltctrl=ControlPlots)

