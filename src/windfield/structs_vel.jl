# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# The different wind velocity types (vel_mode).

"""
    VelModel

Abstract type representing a velocity model in the wind field module.
Subtypes of `VelModel` implement specific velocity field representations or models.

See also: 
- [Defining the wind velocity model](@ref) for more details.
"""
abstract type VelModel end

"""
    Velocity_Constant <: VelModel

A velocity model representing a constant wind velocity field. This struct is used as a type marker 
to indicate that the wind velocity does not vary in space or time.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_Constant <: VelModel end

"""
    Velocity_Constant_wErrorCov

A velocity model representing a constant wind field with associated 
error covariance. This struct is a subtype of `VelModel` and is used 
to model wind velocity with an constant value and an error covariance 
for uncertainty quantification.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_Constant_wErrorCov <: VelModel end

"""
    Velocity_EnKF_InterpTurbine

A velocity model type representing an interpolated turbine velocity field using the 
Ensemble Kalman Filter (EnKF) approach.

# Description
This struct is used within the wind field modeling framework to represent the velocity at a turbine location, 
where the velocity is estimated or interpolated using EnKF-based techniques.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_EnKF_InterpTurbine <: VelModel end
struct Velocity_I_and_I <: VelModel end
struct Velocity_Interpolation <: VelModel end
struct Velocity_Interpolation_wErrorCov <: VelModel end
struct Velocity_InterpTurbine <: VelModel end
struct Velocity_InterpTurbine_wErrorCov <: VelModel end
struct Velocity_RW_with_Mean <: VelModel end
struct Velocity_ZOH_wErrorCov <: VelModel end