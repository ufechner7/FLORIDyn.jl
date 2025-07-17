# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# the different wind velocity types (vel_mode)

struct Velocity_Constant end
struct Velocity_Constant_wErrorCov end
struct Velocity_EnKF_InterpTurbine end
struct Velocity_I_and_I end
struct Velocity_Interpolation end
struct Velocity_Interpolation_wErrorCov end
struct Velocity_InterpTurbine end
struct Velocity_InterpTurbine_wErrorCov end
struct Velocity_RW_with_Mean end
struct Velocity_ZOH_wErrorCov end