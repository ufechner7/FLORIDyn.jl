# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Statistics, LinearAlgebra  # For mean, but optional
t_end = 3780.0  # total simulation time [s]

function calc_wind(time)
    local wind
    t_skip = 0.0
    low_wind = 6.0
    high_wind = 8.2
    t1 = t_skip + T_START  # Time to start increased wind speed
    t2 = t_skip + T_END    # Time to stop  increased wind speed
    if time < t1
        wind = low_wind
    elseif time < t2
        wind = high_wind
    else
        wind = low_wind
    end
    return wind
end

"""
    compute_c(u_entry::Vector{Float64}, u_meas::Vector{Float64}, Δt::Float64, x::Float64) -> Float64

Computes advection speed c = x / τ using cross-correlation lag τ.
u_entry: wind speed time series at entry (x=0)
u_meas: wind speed time series at position x (same length and Δt)
Δt: time step [s]
x: distance [m]

Assumes data aligned in time, computes max correlation lag τ.
"""
function compute_c(u_entry::Vector{Float64}, u_meas::Vector{Float64}, Δt::Float64, x::Float64)
    N = length(u_entry)
    @assert length(u_meas) == N "Time series must have equal length"
    
    # Remove mean for better correlation
    u_entry_demean = u_entry .- mean(u_entry)
    u_meas_demean = u_meas .- mean(u_meas)
    
    # Compute cross-correlation lags (negative lags: entry leads meas)
    max_lag = div(N, 4)  # Limit search range for efficiency
    lags = -max_lag:max_lag
    corrs = zeros(Float64, length(lags))
    
    for (i, lag) in enumerate(lags)
        if lag >= 0
            idx1 = 1:N-lag
            idx2 = 1+lag:N
        else
            idx1 = 1-lag:N
            idx2 = 1:N+lag
        end
        corrs[i] = dot(u_entry_demean[idx1], u_meas_demean[idx2])
    end
    
    # Find lag τ maximizing correlation (negative τ means entry arrives first)
    τ_idx = argmax(corrs)
    τ_steps = lags[τ_idx]
    τ_sec = τ_steps * Δt
    
    c = x / τ_sec  # Should be positive if τ_sec > 0
    return c, τ_sec, τ_steps
end

## Example usage
time_step = 4.0  # s
time_vec = 0:time_step:t_end
u0 = calc_wind.(time_vec)  # Simulated entry speeds
c_true = mean(u0)     # m/s
time_step = 4.0          # s
x = 400.0         # m
delay_steps = round(Int, x / c_true / time_step)  # Number of time steps for delay

# Create delayed signal: u_x[i] = u0[i - delay_steps]
u_x = zeros(Float64, length(u0))
for i in 1:length(u0)
    src_idx = max(1, i - delay_steps)
    u_x[i] = u0[src_idx]
end

c_est, τ, τ_steps = compute_c(u0, u_x, time_step, x)
println("Estimated c: $c_est m/s, true: $c_true m/s, τ: $τ s, delay steps: $delay_steps, estimated steps: $τ_steps")

