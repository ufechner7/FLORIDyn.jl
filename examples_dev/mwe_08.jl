function interpolate_scaling(time::Float64, t1::Float64, t2::Float64, scaling::Vector{Float64})
    scaling_begin = scaling[1]
    scaling_mid = scaling[2]
    scaling_end = scaling[3]
    # Quadratic interpolation that satisfies:
    # time=t1 -> scaling_begin, time=(t1+t2)/2 -> scaling_mid, time=t2 -> scaling_end
    s = clamp((time - t1) / (t2 - t1), 0.0, 1.0)
    L0 = 2 * (s - 0.5) * (s - 1.0)    # basis for value at s=0
    L1 = 4 * s * (1.0 - s)            # basis for value at s=0.5
    L2 = 2 * s * (s - 0.5)            # basis for value at s=1
    scaling = scaling_begin * L0 + scaling_mid * L1 + scaling_end * L2
    scaling = max(scaling_begin, min(scaling_end, scaling))  # clamp scaling to [scaling_begin, scaling_end]
end

# Example usage
t1 = 0.0
t2 = 10.0
scaling = [1.0, 1.2, 2.0]
for time in 0.0:1.0:10.0
    s = interpolate_scaling(time, t1, t2, scaling)
    println("Time: $time, Scaling: $s")
end