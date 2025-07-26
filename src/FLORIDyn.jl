# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

module FLORIDyn

using Interpolations, LinearAlgebra, Random, YAML, StructMapping, Parameters, CSV, DataFrames, DelimitedFiles, JLD2
using Statistics, StaticArrays

export setup, str2type, Settings, getTurbineData, initSimulation

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
export WindVelType, WindVelMatrix, WindFarm

export Direction_All, Direction_Influence, Direction_None, Direction_wGaspariAndCohn
export Velocity_Influence, Velocity_None, Velocity_wGaspariAndCohn
export TI_Influence, TI_None, TI_wGaspariAndCohn
export IterateOPs_average, IterateOPs_basic, IterateOPs_buffer, IterateOPs_maximum, IterateOPs_weighted
export Yaw_Constant, Yaw_InterpTurbine, Yaw_SOWFA

export getWindDirT, getWindDirT_EnKF
export getWindShearT
export getWindTiT, getWindTiT_EnKF
export getWindSpeedT, getWindSpeedT_EnKF
export getDataDir, getDataTI, getDataVel
export correctDir!

export discretizeRotor, calcCt, States
export turbineArrayProperties
export prepareSimulation, importSOWFAFile, centerline, angSOWFA2world, initSimulation
export runFLORIS, init_states, getUadv
export runFLORIDyn, iterateOPs!, getVars, setUpTmpWFAndRun, interpolateOPs, perturbationOfTheWF!, findTurbineGroups

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

# marker structs
include("windfield/structs_dir.jl")
include("windfield/structs_shear.jl")
include("windfield/structs_turb.jl")
include("windfield/structs_vel.jl")
include("correction/structs_dir.jl")
include("correction/structs_vel.jl")
include("correction/structs_turb.jl")
include("floridyn_cl/structs.jl")
include("controller/structs_controller.jl")

"""
    Settings

A mutable struct that holds configuration parameters for the FLORIDyn simulation environment.

# Fields
- `vel_mode::VelModel`: See: [VelModel](@ref)
- `dir_mode::DirModel`: See: [DirModel](@ref)
- `turb_mode`
- `shear_mode`
- `cor_dir_mode`
- `cor_vel_mode`
- `cor_turb_mode`
- `iterate_mode`
- `control_mode`
"""
mutable struct Settings
    vel_mode::VelModel
    dir_mode::DirModel
    turb_mode
    shear_mode
    cor_dir_mode
    cor_vel_mode
    cor_turb_mode
    iterate_mode
    control_mode
end

"""
    WindShear

A struct representing the wind shear profile. This type is used to model the variation of wind speed with height, 
which is important in atmospheric and wind energy simulations.

# Fields
- z0::Float64: Reference height (not used in the [`getWindShearT`](@ref))
- alpha::Float64: WindShear coefficient
"""
struct WindShear
    alpha::Float64
    z0::Float64
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
- `Data`::Matrix{Float64}:    Columns [time, phi] or [time, `phi_T0`, `phi_T1`, ... `phi_Tn`]
- `CholSig`::Matrix{Float64}: Cholesky factor of covariance matrix (nT x nT)
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

"""
    WindFarm

A mutable struct representing a wind farm. Fields can be specified using keyword arguments.

# Fields
- nT::Int64: Number of turbines
- nOP::Int64: Number of operating points
- States_WF::Matrix{Float64}: States of the wind farm
- States_OP::Matrix{Float64}: States of the operating points
- States_T::Matrix{Float64}: States of the turbines
- posBase::Matrix{Float64}: Base positions of the turbines
- posNac::Matrix{Float64}: Positions of the nacelles
- D::Vector{Float64}: Diameters of the turbines
- StartI::Matrix{Int}: Start indices for each turbine
- intOPs::Vector{Matrix{Float64}}: Interpolated operating points
- Weight::Vector{Vector{Float64}}: Weights for the operating points
- dep::Vector{Vector{Int}}: Dependencies between turbines
- red_arr::Matrix{Float64}: Reduced array for each turbine
- Names_T::Vector{String}: Names of the states of the turbines
- Names_WF::Vector{String}: Names of the states of the wind farm
- Names_OP::Vector{String}: Names of coordinates the operating points   
"""
@kwdef mutable struct WindFarm
    nT::Int64 = 0                                               # Number of turbines
    nOP::Int64 = 0                                              # Number of operating points
    States_WF::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)   # States of the wind farm
    States_OP::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)   # States of the operating points
    States_T::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)    # States of the turbines
    posBase::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)     # Base positions of the turbines
    posNac::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)      # Positions of the nacelles
    D::Vector{Float64} = Vector{Float64}(undef, 0)              # Diameters of the turbines
    StartI::Matrix{Int} = Matrix{Int}(undef, 0, 0)              # Start indices for each turbine
    intOPs::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}() # Interpolated operating points
    Weight::Vector{Vector{Float64}} = Vector{Vector{Float64}}() # Weights
    dep::Vector{Vector{Int}} = Vector{Vector{Int}}()            # Dependencies between turbines
    red_arr::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)     # Reduced array for each turbine
    Names_T::Vector{String} = Vector{String}(undef, 0)          # Names of the states of the turbines
    Names_WF::Vector{String} = Vector{String}(undef, 0)         # Names of the states of the wind farm
    Names_OP::Vector{String} = Vector{String}(undef, 0)         # Names of the states of the operating points
end

include("settings.jl")
include("init_turbines.jl")

# functions for calculating the wind field
include("windfield/windfield_direction.jl")
include("windfield/windfield_shear.jl")
include("windfield/windfield_turbulence.jl")
include("windfield/windfield_velocity.jl")

include("floris/discretization.jl")
include("floris/gaussian.jl")
include("floridyn_cl/floridyn_cl.jl")

include("correction/direction.jl")
include("correction/velocity.jl")
include("correction/turbulence.jl")

include("floridyn_cl/prepare_simulation.jl")
include("floridyn_cl/iterate.jl")

include("controller/controller.jl")

end
