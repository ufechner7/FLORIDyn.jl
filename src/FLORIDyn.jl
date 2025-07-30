# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
$(DocStringExtensions.README)
"""
module FLORIDyn

# using ControlPlots

using PrecompileTools: @setup_workload, @compile_workload
using LaTeXStrings
import DocStringExtensions

using Interpolations, LinearAlgebra, Random, YAML, StructMapping, Parameters, CSV, DataFrames, DelimitedFiles, JLD2
using Statistics, StaticArrays, Pkg

export setup, Settings, Vis, getTurbineData, initSimulation, TurbineArray

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
export Floris, FloriDyn, Wind, Sim, Con, IterateOPsBuffers

export Direction_All, Direction_Influence, Direction_None
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
export prepareSimulation, importSOWFAFile, centerline, angSOWFA2world, initSimulation
export runFLORIS, init_states, getUadv
export runFLORIDyn, iterateOPs!, getVars, setUpTmpWFAndRun, interpolateOPs, perturbationOfTheWF!, findTurbineGroups
export getMeasurements, calcFlowField, plotFlowField, plotMeasurements, install_examples

# global variables
RNG::AbstractRNG = Random.default_rng()
function set_rng(rng)
    global RNG
    RNG = rng
end

function str2type(name)
    typename = Symbol(name)
    t = getfield(FLORIDyn, typename)
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
include("visualisation/calc_flowfield.jl")
include("visualisation/plot_flowfield.jl")

"""
    copy_model_settings()

Copy model configuration files and data directories to the local data directory.

This function copies essential model configuration files and simulation data from the package's 
data directory to a local `data/` directory in the current working directory. 

# Files and Directories Copied

## Configuration File
- `2021_9T_Data.yaml`: Main wind farm configuration file containing turbine layout, simulation parameters, 
  and model settings

## Data Directory  
- `2021_9T_Data/`: Complete SOWFA simulation data directory containing:
  - `SOWFA_bladePitch.csv`: Blade pitch angle time series
  - `SOWFA_generatorPower.csv`: Generator power output data
  - `SOWFA_generatorTorque.csv`: Generator torque measurements
  - `SOWFA_nacelleYaw.csv`: Nacelle yaw angle data
  - `SOWFA_rotorSpeedFiltered.csv`: Filtered rotor speed measurements
  - `U.csv`, `WindVel.csv`: Wind velocity data
  - `WindDir.csv`, `WindDirConstant.csv`: Wind direction measurements
  - `WindTI.csv`, `WindTIConstant.csv`: Turbulence intensity data
  - Additional covariance and profile files

# Automatic Operations
The function automatically:
- Creates the `data/` directory if it doesn't exist
- Copies the main YAML configuration file using [`copy_files`](@ref)
- Recursively copies the entire `2021_9T_Data/` subdirectory with all CSV files
- Sets proper file permissions (0o774) on all copied files

This function is called as part of [`install_examples`](@ref) to set up a complete 
working environment with all necessary configuration files and simulation data.

See also: [`copy_files`](@ref), [`install_examples`](@ref)
"""
function copy_model_settings()
    files = ["2021_9T_Data.yaml"]
    dst_path = abspath(joinpath(pwd(), "data"))
    
    # Copy main configuration file
    copy_files("data", files)
    
    # Copy the 2021_9T_Data directory and all its contents
    src_data_dir = joinpath(dirname(pathof(@__MODULE__)), "..", "data", "2021_9T_Data")
    dst_data_dir = joinpath(pwd(), "data", "2021_9T_Data")
    
    if isdir(src_data_dir)
        cp(src_data_dir, dst_data_dir, force=true)
        # Set permissions for all copied files in the directory
        for (root, dirs, files_in_dir) in walkdir(dst_data_dir)
            for file in files_in_dir
                chmod(joinpath(root, file), 0o774)
            end
        end
    end
    
    println("Copied $(length(files)) files and 1 directories to $(dst_path) !")
end

"""
    copy_bin()

Copy the script run_julia to the folder "bin"
(it will be created if it doesn't exist).
"""
function copy_bin()
    PATH = "bin"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    cp(joinpath(src_path, "run_julia"), joinpath(PATH, "run_julia"), force=true)
    chmod(joinpath(PATH, "run_julia"), 0o774)
end

"""
    install_examples(add_packages=true)

Install example files, executables, and data files for the FLORIDyn.jl package.

This function sets up a complete working environment by copying:
- Example Julia scripts from the package's examples directory (via [`copy_examples`](@ref))
- Executable scripts (like `run_julia`) to a local `bin` directory (via [`copy_bin`](@ref))
- Model configuration files and data to a local `data` directory (via [`copy_model_settings`](@ref))

# Arguments
- `add_packages::Bool=true`: Whether to automatically install additional required packages 
  ("LaTeXStrings", "Timers") that are commonly used in the examples

# Example
```julia
# Install examples with automatic package installation
install_examples()

# Install examples without installing additional packages
install_examples(false)
```

After running this function, you can:
- Run example scripts from the created `examples/` directory
- Execute `./bin/run_julia` to start Julia with the project environment
- Access model data files in the `data/` directory

See also: [`copy_examples`](@ref), [`copy_bin`](@ref), [`copy_model_settings`](@ref)
"""
function install_examples(add_packages=true)
    copy_examples()
    copy_bin()
    copy_model_settings()
    if add_packages
        Pkg.add(["LaTeXStrings", "Timers"])
    end
end

"""
    copy_examples()

Copy all example scripts to the folder "examples"
(it will be created if it doesn't exist).
"""
function copy_examples()
    PATH = "examples"
    if ! isdir(PATH) 
        mkdir(PATH)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", PATH)
    copy_files("examples", readdir(src_path))
end

"""
    copy_files(relpath, files)

Copy a list of files from the package directory to a local directory.

This utility function copies files from the package's source directory structure to 
a corresponding local directory. It automatically handles directory creation and 
sets appropriate file permissions.

# Arguments
- `relpath::String`: The relative path (directory name) within both the package and 
  local filesystem where files should be copied
- `files::Vector{String}`: A vector of filenames to copy from the source to destination

# Details
The function:
- Creates the destination directory (`relpath`) if it doesn't exist
- Copies each file from `<package_dir>/<relpath>/<file>` to `./<relpath>/<file>`
- Sets executable permissions (0o774) on all copied files
- Overwrites existing files (uses `force=true`)

# Returns
- `Vector{String}`: The list of files that were copied (same as input `files`)

# Example
```julia
# Copy specific example files to local examples directory
copy_files("examples", ["main.jl", "menu.jl"])

# Copy data files to local data directory  
copy_files("data", ["config.yaml", "turbine_data.csv"])
```

This function is used internally by [`copy_examples`](@ref), [`copy_model_settings`](@ref), 
and other file copying utilities.

See also: [`copy_examples`](@ref), [`copy_model_settings`](@ref), [`copy_bin`](@ref)
"""
function copy_files(relpath, files)
    if ! isdir(relpath) 
        mkdir(relpath)
    end
    src_path = joinpath(dirname(pathof(@__MODULE__)), "..", relpath)
    for file in files
        cp(joinpath(src_path, file), joinpath(relpath, file), force=true)
        chmod(joinpath(relpath, file), 0o774)
    end
    files
end

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster.
    path = dirname(pathof(@__MODULE__))
    path = (joinpath(path, "..", "data"))
    vis = Vis(online=false)
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        settings_file = joinpath(path, "2021_9T_Data.yaml")
        wind, sim, con, floris, floridyn, ta = setup(settings_file)
        set = Settings(wind, sim, con)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        wf = initSimulation(wf, sim)
        runFLORIDyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)
    end

end
end
