# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# The different velocity correction types.

"""
    VelCorrection

An abstract type representing a velocity correction model for wind field adjustments.
Subtypes of `VelCorrection` should implement specific correction methods for wind velocity.
"""
abstract type VelCorrection end

"""
    Velocity_Influence <: VelCorrection

A marker struct used to represent velocity correction based on influence modeling.
"""
struct Velocity_Influence <: VelCorrection end

"""
    Velocity_None <: VelCorrection

A marker struct used to indicate that no velocity corrections should be applied.
"""
struct Velocity_None <: VelCorrection end

"""
    Velocity_wGaspariAndCohn <: VelCorrection

A marker struct used to represent velocity correction using the Gaspari and Cohn localization method.
This correction method is commonly used in ensemble data assimilation for spatial localization.
"""
struct Velocity_wGaspariAndCohn <: VelCorrection end
