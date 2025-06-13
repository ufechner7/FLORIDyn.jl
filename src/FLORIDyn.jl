# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

module FLORIDyn

using Interpolations, LinearAlgebra

export Direction_Constant, Direction_Constant_wErrCov, Direction_EnKF_InterpTurbine, Direction_Interpolation
export Direction_Interpolation_wErrorCov, Direction_InterpTurbine, Direction_InterpTurbine_wErrorCov
export Direction_RW_with_Mean

export getWindDirT, getWindDirT_EnKF

# the different wind direction types (dir_mode)
struct Direction_Constant end
struct Direction_Constant_wErrCov end
struct Direction_EnKF_InterpTurbine end
struct Direction_Interpolation end
struct Direction_Interpolation_wErrorCov end
struct Direction_InterpTurbine end
struct Direction_InterpTurbine_wErrorCov end
struct Direction_RW_with_Mean end

include("windfield.jl")

end
