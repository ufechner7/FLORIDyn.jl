# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Helper function for piecewise cubic Hermite spline interpolation

"""
    interpolate_hermite_spline(s::Float64, correction::Vector, s_positions::Vector) -> Float64

Perform piecewise cubic Hermite spline interpolation across non-equidistant control points.

This function uses cubic Hermite spline interpolation between control points 
located at specified positions along s ∈ [0, 1]. The method provides smooth C1-continuous 
transitions while respecting the control point values.

# Arguments
- `s::Float64`: Normalized parameter in [0, 1] representing position along the curve
- `correction::Vector`: Vector containing control point values
- `s_positions::Vector`: Normalized positions [0, 1] of each control point

# Returns
- `Float64`: Interpolated value at position s
"""
function interpolate_hermite_spline(s::Float64, correction::Vector, s_positions::Vector)
    n_points = length(correction)
    @assert n_points >= 2 "Need at least 2 control points for interpolation"
    @assert length(s_positions) == n_points "Number of positions must match number of corrections"
    @assert s_positions[1] == 0.0 && s_positions[end] == 1.0 "Positions must span [0, 1]"
    
    # Handle edge cases
    if n_points == 2
        # Linear interpolation for 2 points
        return correction[1] + s * (correction[2] - correction[1])
    end
    
    # Find which segment we're in
    segment_idx = 1
    for i in 1:(n_points-1)
        if s >= s_positions[i] && s <= s_positions[i+1]
            segment_idx = i
            break
        end
    end
    if s >= s_positions[end]
        segment_idx = n_points - 1
    end
    
    # Calculate tangents for all control points using finite differences
    tangents = zeros(n_points)
    
    # First point: forward difference
    tangents[1] = (correction[2] - correction[1]) / (s_positions[2] - s_positions[1])
    
    # Interior points: central differences (accounting for non-uniform spacing)
    for i in 2:(n_points-1)
        h_left = s_positions[i] - s_positions[i-1]
        h_right = s_positions[i+1] - s_positions[i]
        # Weighted average based on segment widths
        tangents[i] = ((correction[i] - correction[i-1]) / h_left + 
                       (correction[i+1] - correction[i]) / h_right) / 2.0
    end
    
    # Last point: backward difference
    tangents[n_points] = (correction[n_points] - correction[n_points-1]) / 
                         (s_positions[n_points] - s_positions[n_points-1])
    
    # Local parameter t within the segment [0, 1]
    s_start = s_positions[segment_idx]
    s_end = s_positions[segment_idx + 1]
    segment_width = s_end - s_start
    t = clamp((s - s_start) / segment_width, 0.0, 1.0)
    
    # Cubic Hermite basis functions
    h00 = 2*t^3 - 3*t^2 + 1
    h10 = t^3 - 2*t^2 + t
    h01 = -2*t^3 + 3*t^2
    h11 = t^3 - t^2
    
    # Interpolate using values and tangents at segment endpoints
    p0 = correction[segment_idx]
    p1 = correction[segment_idx + 1]
    m0 = tangents[segment_idx]
    m1 = tangents[segment_idx + 1]
    
    correction_result = h00*p0 + h10*segment_width*m0 + h01*p1 + h11*segment_width*m1
    
    return correction_result
end

"""
    interpolate_hermite_spline(s::Float64, correction::Vector) -> Float64

Perform piecewise cubic Hermite spline interpolation across control points.

This function uses cubic Hermite spline interpolation between control points 
evenly spaced along s ∈ [0, 1]. The method provides smooth C1-continuous transitions while
respecting the control point values. The number of control points is determined
from the length of the `correction` vector.

# Arguments
- `s::Float64`: Normalized parameter in [0, 1] representing position along the curve
- `correction::Vector`: Vector containing control point values. For n control points,
  they are located at s = 0, 1/(n-1), 2/(n-1), ..., 1.0

# Returns
- `Float64`: Interpolated value at position s
"""
function interpolate_hermite_spline(s::Float64, correction::Vector)
    n_points = length(correction)
    @assert n_points >= 2 "Need at least 2 control points for interpolation"
    
    # Handle edge cases
    if n_points == 2
        # Linear interpolation for 2 points
        return correction[1] + s * (correction[2] - correction[1])
    end
    
    # Number of segments = number of control points - 1
    n_segments = n_points - 1
    segment_width = 1.0 / n_segments
    
    # Calculate tangents for all control points using finite differences
    tangents = zeros(n_points)
    
    # First point: forward difference
    tangents[1] = (correction[2] - correction[1]) / segment_width
    
    # Interior points: central differences
    for i in 2:(n_points-1)
        tangents[i] = (correction[i+1] - correction[i-1]) / (2 * segment_width)
    end
    
    # Last point: backward difference
    tangents[n_points] = (correction[n_points] - correction[n_points-1]) / segment_width
    
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
    p0 = correction[segment_idx]
    p1 = correction[segment_idx + 1]
    m0 = tangents[segment_idx]
    m1 = tangents[segment_idx + 1]
    
    correction_result = h00*p0 + h10*segment_width*m0 + h01*p1 + h11*segment_width*m1
    
    return correction_result
end
