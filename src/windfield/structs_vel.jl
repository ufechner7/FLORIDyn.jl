# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# the different wind velocity types (vel_mode)

"""
    VelModel

Abstract type representing a velocity model in the wind field module.
Subtypes of `VelModel` implement specific velocity field representations or models.
"""
abstract type VelModel end

struct Velocity_Constant <: VelModel end
struct Velocity_Constant_wErrorCov <: VelModel end
struct Velocity_EnKF_InterpTurbine <: VelModel end
struct Velocity_I_and_I <: VelModel end
struct Velocity_Interpolation <: VelModel end
struct Velocity_Interpolation_wErrorCov <: VelModel end
struct Velocity_InterpTurbine <: VelModel end
struct Velocity_InterpTurbine_wErrorCov <: VelModel end
struct Velocity_RW_with_Mean <: VelModel end
struct Velocity_ZOH_wErrorCov <: VelModel end