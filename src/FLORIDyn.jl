# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

module FLORIDyn

using Interpolations, LinearAlgebra, Random, YAML, StructMapping, Parameters

export setup, str2type, Settings

export Direction_Constant, Direction_Constant_wErrorCov, Direction_EnKF_InterpTurbine, Direction_Interpolation
export Direction_Interpolation_wErrorCov, Direction_InterpTurbine, Direction_InterpTurbine_wErrorCov
export Direction_RW_with_Mean
export Shear_Interpolation, Shear_PowerLaw, WindShear
export TI_Constant, TI_EnKF_InterpTurbine, TI_Interpolation, TI_InterpTurbine
export Velocity_Constant, Velocity_Constant_wErrorCov, Velocity_EnKF_InterpTurbine
export Velocity_I_and_I, Velocity_Interpolation, Velocity_Interpolation_wErrorCov
export Velocity_InterpTurbine, Velocity_InterpTurbine_wErrorCov, Velocity_RW_with_Mean
export Velocity_ZOH_wErrorCov

export WindDirType, WindDirMatrix, WindDirTriple
export WindVelType, WindVelMatrix

export Direction_All, Direction_Influence, Direction_None, Direction_wGaspariAndCohn
export Velocity_Influence, Velocity_None, Velocity_wGaspariAndCohn
export TI_Influence, TI_None, TI_wGaspariAndCohn
export IterateOPs_average, IterateOPs_basic, IterateOPs_buffer, IterateOPs_maximum, IterateOPs_weighted

export getWindDirT, getWindDirT_EnKF
export getWindShearT
export getWindTiT, getWindTiT_EnKF
export getWindSpeedT, getWindSpeedT_EnKF

export discretizeRotor, CalcCt

# global variables
RNG::AbstractRNG = Random.default_rng()
function set_rng(rng)
    global RNG
    RNG = rng
end

function str2type(name)
    typename = Symbol(name)
    t = getfield(Main, typename)
    instance = t()
end

struct Settings
    vel_mode
    dir_mode
    turb_mode
    shear_mode
    cor_dir_mode
    cor_vel_mode
    cor_turb_mode
    iterate_mode
end

# marker structs
include("windfield/structs_dir.jl")
include("windfield/structs_shear.jl")
include("windfield/structs_turb.jl")
include("windfield/structs_vel.jl")
include("correction/structs_dir.jl")
include("correction/structs_vel.jl")
include("correction/structs_turb.jl")
include("floridyn_cl/structs.jl")

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

"""
     WindDirType

# Fields
- Data::Float64: wind direction value
- CholSig::Matrix{Float64}: Cholesky factor of covariance matrix (nT x nT)
"""
struct WindDirType
    Data::Float64
    CholSig::Matrix{Float64}
end

"""
     WindVelType

# Fields
- Data::Float64: wind speed
- CholSig::Matrix{Float64}: Cholesky factor of covariance matrix (nT x nT)
"""
struct WindVelType
    Data::Float64
    CholSig::Matrix{Float64}
end

"""
    struct WindDirMatrix

# Fields
- Data::Matrix{Float64}:    Nx2 matrix: column 1 = time, column 2 = phi
- CholSig::Matrix{Float64}: Cholesky factor of covariance matrix (nT x nT)
"""
struct WindDirMatrix
    Data::Matrix{Float64}      # Nx2 matrix: column 1 = time, column 2 = phi
    CholSig::Matrix{Float64}   # Cholesky factor of covariance matrix (nT x nT)
end

"""
    struct WindVelMatrix

# Fields
- Data::Matrix{Float64}:    Nx2 matrix: column 1 = time, column 2 = wind speed
- CholSig::Matrix{Float64}: Cholesky factor of covariance matrix (nT x nT)
"""
struct WindVelMatrix
    Data::Matrix{Float64}
    CholSig::Matrix{Float64}
end

"""
    WindDirTriple

A structure representing a wind direction triple. 

# Fields
- Init::Vector{Float64}:    Mean direction (vector or scalar)
- CholSig::Matrix{Float64}: Cholesky factor of covariance matrix (nT x nT)
- MeanPull::Float64:        Scalar mean reversion factor
"""
struct WindDirTriple
    Init::Vector{Float64}      # Mean direction (vector or scalar)
    CholSig::Matrix{Float64}   # Cholesky factor of covariance matrix (nT x nT)
    MeanPull::Float64          # Scalar mean reversion factor
end

include("settings.jl")
# functions for calculating the wind field
include("windfield/windfield_direction.jl")
include("windfield/windfield_shear.jl")
include("windfield/windfield_turbulence.jl")
include("windfield/windfield_velocity.jl")

include("floris/discretization.jl")
include("floris/gaussian.jl")

include("correction/direction.jl")
include("correction/velocity.jl")

end
