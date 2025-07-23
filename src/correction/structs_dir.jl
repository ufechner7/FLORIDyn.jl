# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# The different direction correction types.

"""
    DirCorrection

An abstract type representing a direction correction model for wind field adjustments.
Subtypes of `DirCorrection` should implement specific correction methods for wind direction.

See: [Markers defining the direction correction](@ref) for more details.
"""
abstract type DirCorrection end

"""
    Direction_All <: DirCorrection

A marker struct used to indicate that all direction corrections should be applied.
"""
struct Direction_All <: DirCorrection end

"""
    Direction_Influence <: DirCorrection

A marker struct used to represent direction correction based on influence modeling.
"""
struct Direction_Influence <: DirCorrection end

"""
    Direction_None <: DirCorrection

A marker struct used to indicate that no direction corrections should be applied.
"""
struct Direction_None <: DirCorrection end

"""
    Direction_wGaspariAndCohn <: DirCorrection

A marker struct used to represent direction correction using the Gaspari and Cohn localization method.
This correction method is commonly used in ensemble data assimilation for spatial localization.
"""
struct Direction_wGaspariAndCohn <: DirCorrection end
