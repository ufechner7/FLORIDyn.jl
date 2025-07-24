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
"""
struct TI_Influence <: TurbulenceCorrection end

"""
    TI_None <: TurbulenceCorrection

A marker struct used to indicate that no turbulence intensity corrections should be applied.
"""
struct TI_None <: TurbulenceCorrection end

"""
    TI_wGaspariAndCohn <: TurbulenceCorrection

A marker struct used to represent turbulence intensity correction using the Gaspari and Cohn localization method.
This correction method is commonly used in ensemble data assimilation for spatial localization.
"""
struct TI_wGaspariAndCohn <: TurbulenceCorrection end