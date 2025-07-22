# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# the different wind direction types (dir_mode)

"""
    DirModel

An abstract type representing a directional wind field model. 
Subtypes of `DirModel` should implement specific models for wind directionality.

See: [Markers defining the wind direction model](@ref) for more details.
"""
abstract type DirModel end

"""
    Direction_Constant <: DirModel

A marker struct used to represent a constant wind direction.

## Example:
```julia
dir_mode = Direction_constant()
phi = getWindDirT(dir_mode, 270, [1,2,3], nothing)
```
"""
struct Direction_Constant <: DirModel end
"""
    Direction_Constant_wErrorCov <: DirModel

A marker struct used to indicate a wind direction that is constant with associated error covariance.
"""
struct Direction_Constant_wErrorCov <: DirModel end
"""
    Direction_EnKF_InterpTurbine <: DirModel

A marker struct used to indicate the use of direction-aware Ensemble Kalman Filter (EnKF) interpolation for turbine modeling.
"""
struct Direction_EnKF_InterpTurbine <: DirModel end
"""
    Direction_Interpolation <: DirModel

A marker struct used to represent direction interpolation functionality within the FLORIDyn framework.
"""
struct Direction_Interpolation <: DirModel end
"""
    Direction_Interpolation_wErrorCov <: DirModel

A marker struct representing a direction interpolation method with associated error covariance.
"""
struct Direction_Interpolation_wErrorCov <: DirModel end
"""
    Direction_InterpTurbine <: DirModel

A marker struct used to indicate direction interpolation for turbines.
"""
struct Direction_InterpTurbine <: DirModel end
"""
    Direction_InterpTurbine_wErrorCov <: DirModel

A marker struct used to indicate the use of direction interpolation for turbines with associated error covariance.
"""
struct Direction_InterpTurbine_wErrorCov <: DirModel end
"""
    Direction_RW_with_Mean <: DirModel

A marker struct used to indicate the use of a random walk direction model with a mean component.
"""
struct Direction_RW_with_Mean <: DirModel end
