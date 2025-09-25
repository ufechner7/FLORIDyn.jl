# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# The different turbulence correction types.

"""
    TurbulenceCorrection

An abstract type representing a turbulence correction model for wind field adjustments.
Subtypes of `TurbulenceCorrection` should implement specific correction methods for turbulence intensity.

See: [Defining the turbulence correction](@ref) for more details.
"""
abstract type TurbulenceCorrection end

"""
    TI_Influence <: TurbulenceCorrection

A marker struct used to represent turbulence intensity correction based on influence modeling.

# WARNING:
This correction method is not properly tested. Use at your own risk!
"""
struct TI_Influence <: TurbulenceCorrection end

"""
    TI_None <: TurbulenceCorrection

A marker struct used to indicate that no turbulence intensity corrections should be applied.
"""
struct TI_None <: TurbulenceCorrection end
