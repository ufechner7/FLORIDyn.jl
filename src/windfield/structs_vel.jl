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
    Velocity_Constant_wErrorCov <: VelModel

A velocity model representing a constant wind field with associated 
error covariance. This struct is a subtype of `VelModel` and is used 
to model wind velocity with an constant value and an error covariance 
for uncertainty quantification.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_Constant_wErrorCov <: VelModel end

"""
    Velocity_EnKF_InterpTurbine <: VelModel

A velocity model type representing an interpolated turbine velocity field using the 
Ensemble Kalman Filter (EnKF) approach.

# Description
This struct is used within the wind field modeling framework to represent the velocity at a turbine location, 
where the velocity is estimated or interpolated using EnKF-based techniques.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_EnKF_InterpTurbine <: VelModel end

"""
    Velocity_I_and_I <: VelModel

A velocity model implementing an interpolation and integration approach for wind velocity estimation.

# Description
This struct represents a velocity model that combines interpolation techniques with integration methods 
to estimate wind velocity fields, typically used for advanced wind field reconstruction scenarios.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_I_and_I <: VelModel end

"""
    Velocity_Interpolation <: VelModel

A velocity model that uses spatial interpolation techniques to estimate wind velocity fields.

# Description
This struct represents a velocity model that employs interpolation methods to determine wind velocities 
at arbitrary spatial locations based on available measurement data or model predictions.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.

# WARNING
This model is **not yet implemented**!
"""
struct Velocity_Interpolation <: VelModel end

"""
    Velocity_Interpolation_wErrorCov <: VelModel

A velocity model that uses spatial interpolation with associated error covariance information.

# Description
This struct represents a velocity model that employs interpolation methods to determine wind velocities 
and includes error covariance matrices for uncertainty quantification and probabilistic analysis.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_Interpolation_wErrorCov <: VelModel end

"""
    Velocity_InterpTurbine <: VelModel

A velocity model for interpolating wind velocities specifically at turbine locations.

# Description
This struct represents a velocity model that focuses on estimating wind velocities at turbine hub heights 
and rotor positions using interpolation techniques from surrounding measurement points or model data.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_InterpTurbine <: VelModel end

"""
    Velocity_InterpTurbine_wErrorCov <: VelModel

A velocity model for interpolating wind velocities at turbine locations with error covariance information.

# Description
This struct represents a velocity model that estimates wind velocities at turbine positions using 
interpolation techniques and includes associated error covariance matrices for uncertainty analysis 
and robust wind farm control applications.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_InterpTurbine_wErrorCov <: VelModel end

"""
    Velocity_RW_with_Mean <: VelModel

A velocity model implementing a random walk process with a mean trend component.

# Description
This struct represents a velocity model that combines a random walk stochastic process with a 
deterministic mean component, typically used for modeling wind velocity evolution over time 
with both predictable trends and random fluctuations.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_RW_with_Mean <: VelModel end

"""
    Velocity_ZOH_wErrorCov <: VelModel

A velocity model using Zero-Order Hold (ZOH) interpolation with error covariance information.

# Description
This struct represents a velocity model that employs zero-order hold interpolation (piecewise constant) 
for wind velocity estimation between measurement points, and includes error covariance matrices for 
uncertainty quantification and statistical analysis.

# See also:
- [`VelModel`](@ref): Abstract supertype for velocity models.
"""
struct Velocity_ZOH_wErrorCov <: VelModel end