# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# the different wind velocity types (vel_mode)

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