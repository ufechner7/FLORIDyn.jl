# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# the different wind turbulence types (turb_mode)

"""
    TurbulenceModel

An abstract type representing a turbulence model for wind field calculations.
Subtypes of `TurbulenceModel` should implement specific models for wind turbulence intensity.

See: [Markers for defining the wind turbulence](@ref) for more details.
"""
abstract type TurbulenceModel end

"""
    TI_Constant <: TurbulenceModel

A marker struct representing a constant turbulence intensity. 
"""
struct TI_Constant <: TurbulenceModel end

"""
    TI_EnKF_InterpTurbine <: TurbulenceModel

A marker struct representing the Turbulence Intensity (TI) Ensemble Kalman Filter (EnKF) interpolation model.
"""
struct TI_EnKF_InterpTurbine <: TurbulenceModel end

"""
    TI_Interpolation <: TurbulenceModel

A marker struct representing the interpolation method for modeling the turbulence.
"""
struct TI_Interpolation <: TurbulenceModel end


"""
    TI_InterpTurbine <: TurbulenceModel

A marker struct representing an interpolated turbine model for turbulence intensity calculations.
"""
struct TI_InterpTurbine <: TurbulenceModel end
