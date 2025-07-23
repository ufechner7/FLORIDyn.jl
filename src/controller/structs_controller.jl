# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# The different controller types.

"""
    ControllerModel

An abstract type representing a controller model for wind turbines.
Subtypes of `ControllerModel` should implement specific control strategies for turbine operation.

See: [Markers defining the controller](@ref) for more details.
"""
abstract type ControllerModel end

"""
    Yaw_Constant <: ControllerModel

A marker struct used to represent a constant yaw control strategy.
In this mode, turbines maintain a fixed yaw angle throughout the simulation.
"""
struct Yaw_Constant <: ControllerModel end

"""
    Yaw_InterpTurbine <: ControllerModel

A marker struct used to indicate yaw control with turbine interpolation.
This mode allows for interpolated yaw angles across different turbine positions.
"""
struct Yaw_InterpTurbine <: ControllerModel end

"""
    Yaw_SOWFA <: ControllerModel

A marker struct used to represent yaw control compatible with SOWFA (Simulator fOr Wind Farm Applications).
This mode is specifically designed for integration with SOWFA simulation data.
"""
struct Yaw_SOWFA <: ControllerModel end