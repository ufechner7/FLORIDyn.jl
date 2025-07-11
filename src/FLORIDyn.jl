# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

module FLORIDyn

using Interpolations, LinearAlgebra, Random

export Direction_Constant, Direction_Constant_wErrorCov, Direction_EnKF_InterpTurbine, Direction_Interpolation
export Direction_Interpolation_wErrorCov, Direction_InterpTurbine, Direction_InterpTurbine_wErrorCov
export Direction_RW_with_Mean
export Shear_Interpolation, Shear_PowerLaw, WindShear
export TI_Constant, TI_EnKF_InterpTurbine, TI_Interpolation, TI_InterpTurbine

export getWindDirT, getWindDirT_EnKF
export getWindShearT
export getWindTiT, getWindTiT_EnKF

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

"""
    Shear_Interpolation

A marker struct used to represent the linear interpolation for wind shear profiles.

# See also
- [`Shear_PowerLaw`](@ref)
"""
struct Shear_Interpolation end

# """
#     Shear_LogLaw

# A type representing the logarithmic law for modeling wind shear profiles.

# # See also
# - [`Shear_Interpolation`](@ref)
# - [`Shear_PowerLaw`](@ref)
# """
# struct Shear_LogLaw end

"""
    Shear_PowerLaw

A marker struct representing the logarithmic law for modeling wind shear profiles.

# See also
- [`Shear_Interpolation`](@ref)
"""
struct Shear_PowerLaw end

"""
    TI_Constant

A marker struct representing a constant turbulence intensity. 
"""
struct TI_Constant end

"""
    TI_EnKF_InterpTurbine

A marker struct representing the Turbulence Intensity (TI) Ensemble Kalman Filter (EnKF) interpolation model.
"""
struct TI_EnKF_InterpTurbine end

"""
    TI_Interpolation

A marker struct representing the interpolation method for modeling the turbulence.
"""
struct TI_Interpolation end
struct TI_InterpTurbine end

function set_rng(rng)
    global RNG
    RNG = rng
end

"""
    WindShear

A struct representing the wind shear profile. This type is used to model the variation of wind speed with height, 
which is important in atmospheric and wind energy simulations.

# Fields
- z0::Float64: Reference height (not used in the [`getWindShearT`](@ref))
- alpha::Float64: Shear coefficient
"""
struct WindShear
    z0::Float64
    alpha::Float64
end

include("windfield_interpolation.jl")
include("windfield_shear.jl")
include("windfield_turbulence.jl")

end
