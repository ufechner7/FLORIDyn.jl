# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# The different velocity correction types.

"""
    VelCorrection

An abstract type representing a velocity correction model for wind field adjustments.
Subtypes of `VelCorrection` should implement specific correction methods for wind velocity.

See: [Defining the velocity correction](@ref) for more details.
"""
abstract type VelCorrection end

"""
    Velocity_Influence <: VelCorrection

Marker struct selecting the influence-based free-stream velocity correction implemented in
`correctVel!(::Velocity_Influence, ...)`. See that function's docstring for algorithm details.

# WARNING:
This correction method is not properly tested. Use at your own risk!
"""
struct Velocity_Influence <: VelCorrection end

"""
    Velocity_None <: VelCorrection

A marker struct used to indicate that no velocity corrections should be applied.
"""
struct Velocity_None <: VelCorrection end

