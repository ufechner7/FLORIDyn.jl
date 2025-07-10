# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

module FLORIDyn

using Interpolations, LinearAlgebra, Random

export Direction_Constant, Direction_Constant_wErrorCov, Direction_EnKF_InterpTurbine, Direction_Interpolation
export Direction_Interpolation_wErrorCov, Direction_InterpTurbine, Direction_InterpTurbine_wErrorCov
export Direction_RW_with_Mean

export getWindDirT, getWindDirT_EnKF
export getWindShearT

# global variables
RNG::AbstractRNG = Random.default_rng()

# the different wind direction types (dir_mode)
"""
    Direction_Constant

A marker struct used to represent a constant wind direction.

## Example:
```julia
dir_mode = Direction_constant()
phi = getWindDirT(dir_mode, 270, [1,2,3], nothing)
```
"""
struct Direction_Constant end
"""
    Direction_Constant_wErrorCov

A marker struct used to indicate a wind direction that is constant with associated error covariance.
"""
struct Direction_Constant_wErrorCov end
"""
    Direction_EnKF_InterpTurbine

A marker struct used to indicate the use of direction-aware Ensemble Kalman Filter (EnKF) interpolation for turbine modeling.
"""
struct Direction_EnKF_InterpTurbine end
"""
    Direction_Interpolation

A marker struct used to represent direction interpolation functionality within the FLORIDyn framework.
"""
struct Direction_Interpolation end
"""
    Direction_Interpolation_wErrorCov

A marker struct representing a direction interpolation method with associated error covariance.
"""
struct Direction_Interpolation_wErrorCov end
"""
    Direction_InterpTurbine

A marker struct used to indicate direction interpolation for turbines.
"""
struct Direction_InterpTurbine end
"""
    Direction_InterpTurbine_wErrorCov

A marker struct used to indicate the use of direction interpolation for turbines with associated error covariance.
"""
struct Direction_InterpTurbine_wErrorCov end
"""
    Direction_RW_with_Mean

A marker struct used to indicate the use of a random walk direction model with a mean component.
"""
struct Direction_RW_with_Mean end

function set_rng(rng)
    global RNG
    RNG = rng
end

include("windfield_interpolation.jl")
include("windfield_shear.jl")

end
