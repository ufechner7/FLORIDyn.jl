# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# The different direction correction types.

"""
    DirCorrection

An abstract type representing a direction correction model for wind field adjustments.
Subtypes of `DirCorrection` should implement specific correction methods for wind direction.

See: [Defining the direction correction](@ref) for more details.
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
