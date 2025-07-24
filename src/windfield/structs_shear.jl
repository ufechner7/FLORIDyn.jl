# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# the different shear types (shear_mode)

"""
    ShearModel

An abstract type representing a wind shear model for vertical wind profiles.
Subtypes of `ShearModel` should implement specific models for wind shear calculations.

See: [Defining the wind shear model](@ref) for more details.
"""
abstract type ShearModel end

"""
    Shear_Interpolation <: ShearModel

A marker struct used to represent the linear interpolation for wind shear profiles.

# See also
- [`Shear_PowerLaw`](@ref)
"""
struct Shear_Interpolation <: ShearModel end

# """
#     Shear_LogLaw

# A type representing the logarithmic law for modeling wind shear profiles.

# # See also
# - [`Shear_Interpolation`](@ref)
# - [`Shear_PowerLaw`](@ref)
# """
# struct Shear_LogLaw end

"""
    Shear_PowerLaw <: ShearModel

A marker struct representing the logarithmic law for modeling wind shear profiles.

# See also
- [`Shear_Interpolation`](@ref)
"""
struct Shear_PowerLaw <: ShearModel end
