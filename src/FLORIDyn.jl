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
function set_rng(rng)
    global RNG
    RNG = rng
end

# marker structs
include("structs_dir.jl")
include("structs_shear.jl")
include("structs_turb.jl")

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

include("windfield_direction.jl")
include("windfield_shear.jl")
include("windfield_turbulence.jl")
include("windfield_velocity.jl")

end
