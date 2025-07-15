# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# the different shear types (shear_mode)
"""
    Shear_Interpolation

A marker struct used to represent the linear interpolation for wind shear profiles.

# See also
- [`Shear_PowerLaw`](@ref)
"""
struct Shear_Interpolation end

# """
#     Shear_LogLaw

# A type representing the logarithmic law for modeling wind shear profiles.

# # See also
# - [`Shear_Interpolation`](@ref)
# - [`Shear_PowerLaw`](@ref)
# """
# struct Shear_LogLaw end

"""
    Shear_PowerLaw

A marker struct representing the logarithmic law for modeling wind shear profiles.

# See also
- [`Shear_Interpolation`](@ref)
"""
struct Shear_PowerLaw end
