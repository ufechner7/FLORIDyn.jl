# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
$(DocStringExtensions.README)
"""
module FLORIDyn

using PrecompileTools: @setup_workload, @compile_workload
import DocStringExtensions
using Interpolations, LinearAlgebra, Random, YAML, StructMapping, Parameters, CSV, DataFrames, DelimitedFiles, JLD2
using Statistics, StaticArrays

export setup, Settings, getTurbineData, initSimulation

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
export getWindTiT
export getWindSpeedT, getWindSpeedT_EnKF
export getDataDir, getDataTI, getDataVel
export correctDir!
export getYaw

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

function str2type(name)
    typename = Symbol(name)
    t = getfield(FLORIDyn, typename)
    instance = t()
end

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

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        settings_file = "data/2021_9T_Data.yaml"
        wind, sim, con, floris, floridyn = setup(settings_file)
        a = Velocity_Constant()
        b = str2type("Velocity_Constant")
        set = Settings(wind, sim, con)
    end
end
end # module FLORIDyn

# --> Without compile_workload:
# Time elapsed: 1.115915423 s
# Time elapsed: 3.398632348 s
# Time elapsed: 3.559121232 s
# Time elapsed: 7.671075958 s
# Time elapsed: 8.222971185 s
#   0.308294 seconds (1.16 M allocations: 190.392 MiB, 8.79% gc time, 66.92% compilation time)

# --> With compile_workload:
# Time elapsed: 1.117815199 s
# Time elapsed: 1.237678098 s
# Time elapsed: 1.387095712 s
# Time elapsed: 5.469984983 s
# Time elapsed: 6.034418235 s
#   0.311222 seconds (1.16 M allocations: 190.379 MiB, 10.51% gc time, 66.15% compilation time)