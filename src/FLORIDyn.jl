# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
$(DocStringExtensions.README)
"""
module FLORIDyn

using PrecompileTools: @setup_workload, @compile_workload
using LaTeXStrings
import DocStringExtensions, LoggingExtras
import Base: show

using Interpolations, LinearAlgebra, Random, YAML, StructMapping, Parameters, CSV, DataFrames, DelimitedFiles, JLD2
using Statistics, StaticArrays, Pkg, DistributedNext, Dates
using REPL.TerminalMenus

export MSR, toMSR, VelReduction, AddedTurbulence, EffWind
export setup, Settings, Vis, getTurbineData, initSimulation, TurbineArray, TurbineData, turbine_group
export set_yaw!, set_induction!

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
export Velocity_Influence, Velocity_None
export TI_Influence, TI_None
export IterateOPs_average, IterateOPs_basic, IterateOPs_buffer, IterateOPs_maximum, IterateOPs_weighted
export Yaw_Constant, Yaw_InterpTurbine, Yaw_SOWFA
export Induction_Constant, Induction_MPC

export getWindDirT, getWindDirT_EnKF
export getWindShearT
export getWindTiT
export getWindSpeedT, getWindSpeedT_EnKF
export getDataDir, getDataTI, getDataVel
export correctDir!, correctTI!, correctVel!
export getYaw, getInduction

export discretizeRotor, calcCt, States
export prepareSimulation, importSOWFAFile, centerline!, angSOWFA2world, initSimulation
export runFLORIS!, init_states, getUadv
export runFLORIDyn, iterateOPs!, setUpTmpWFAndRun!, interpolateOPs!, perturbationOfTheWF!, findTurbineGroups
export getVars!
export getMeasurements, calcFlowField, plotFlowField, plotMeasurements, get_layout, install_examples, calc_rel_power
export run_floridyn, plot_flow_field, plot_measurements, plot_x, plot_rmt, close_all, turbines
export createVideo, createAllVideos, natural_sort_key, cleanup_video_folder
export now_microseconds, now_nanoseconds, precise_now, unique_name, delete_results, find_floridyn_runs, compare_dataframes
export isdelftblue, Measurement, parse_measurements
export FlowField, parse_flow_fields
export UnifiedBuffers, create_unified_buffers
export get_default_project
export select_project
export get_default_msr, set_default_msr, select_measurement

"""
    MSR `VelReduction` `AddedTurbulence` `EffWind`

Enumeration that selects which (scalar) quantity is visualised / stored when plotting
or saving flow field measurements. The acronym stands for Measurement System Representation.
Passed via the `msr` keyword to [`run_floridyn`](@ref) and the plotting helpers
[`plot_flow_field`](@ref), [`plot_measurements`](@ref).

# Elements
- `VelReduction`    (1): Velocity reduction (1 - u / u_ref) downstream of turbines.
- `AddedTurbulence` (2): Added turbulence intensity contributed by wakes (ΔTI component).
- `EffWind`         (3): Effective wind speed at turbine locations (including wake effects).

# Usage
```julia
# Use default (VelReduction)
run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)

# Explicitly request added turbulence visualisation
run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris; msr=AddedTurbulence)

# Convert from a user string (e.g. parsed CLI / YAML value)
msr = toMSR("EffWind")
run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris; msr)
```

See also: [`toMSR`](@ref)
"""
@enum MSR begin
    VelReduction = 1 
    AddedTurbulence = 2 
    EffWind = 3
end

@doc """
VelReduction::MSR
Velocity Reduction.
See also [`MSR`](@ref).
""" VelReduction

@doc """
AddedTurbulence::MSR
Added Turbulence.
See also [`MSR`](@ref).
""" AddedTurbulence

@doc """
EffWind::MSR
Added Turbulence.
See also [`MSR`](@ref).
""" EffWind

"""
    toMSR(s::String)

Converts the input string `s` to a MSR (Measurement System Representation) enumeration.

Supports the following string formats:
- Enum names: `"VelReduction"`, `"AddedTurbulence"`, `"EffWind"`
- Flow field names: `"flow_field_vel_reduction"`, `"flow_field_added_turbulence"`, `"flow_field_eff_wind_speed"`
- Measurement names: `"msr_vel_reduction"`, `"msr_added_turbulence"`, `"msr_eff_wind_speed"`

See also [`MSR`](@ref).
"""
function toMSR(s::String)
    if s == "VelReduction"
        return VelReduction
    elseif s == "AddedTurbulence"
        return AddedTurbulence
    elseif s == "EffWind"
        return EffWind
    # Handle flow field names
    elseif s == "flow_field_vel_reduction"
        return VelReduction
    elseif s == "flow_field_added_turbulence"
        return AddedTurbulence
    elseif s == "flow_field_eff_wind_speed"
        return EffWind
    # Handle measurement names
    elseif s == "msr_vel_reduction"
        return VelReduction
    elseif s == "msr_added_turbulence"
        return AddedTurbulence
    elseif s == "msr_eff_wind_speed"
        return EffWind
    else
        error("Unknown measurement type: $s")
    end
end

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
include("floris/structs_floris.jl")  # Include FLORISBuffers before structs.jl needs it
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
- `parallel::Bool`:  Run plotting in a separate process.
- `threading::Bool`: Enable threading for parallel computation within a single process
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
    induction_mode
    parallel::Bool
    threading::Bool
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

A structure containing the three key parameters for wind direction modeling with mean reversion.

This struct is specifically designed for the `Direction_RW_with_Mean` wind direction mode, 
which implements a random walk with mean reversion model. The model captures both 
stochastic wind direction variability and the tendency for wind to return to prevailing 
climatological directions over time.

# Fields
- `Init::Vector{Float64}`: Target/equilibrium wind directions for each turbine [°].
  These represent the long-term average or expected wind directions 
  that the system tends to revert toward. Each element corresponds to one turbine.
  
- `CholSig::Matrix{Float64}`: Cholesky decomposition of the covariance matrix (nT × nT).
  This matrix encodes both the magnitude of random fluctuations and the spatial 
  correlations between turbines. It allows for:
  - Different noise levels for different turbines (diagonal elements)
  - Cross-correlations between nearby turbines (off-diagonal elements)
  
- `MeanPull::Float64`: Mean reversion strength parameter [0, 1].
  Controls how strongly wind directions are pulled back toward their target values:
  - `0.0`: No mean reversion (pure random walk)
  - `1.0`: Maximum reversion (immediate return to target)
  - Typical values: `0.1` to `0.3` for realistic wind behavior

# Mathematical Model
The wind direction evolves according to:
    φ(t+1) = φ(t) + CholSig × ε + MeanPull × (Init - φ(t))
where `ε` is a vector of independent standard normal random variables.

# Physical Interpretation
- **Meteorological realism**: `Init` represents seasonal/geographic wind patterns
- **Spatial correlation**: `CholSig` models how wind direction changes are correlated across the wind farm
- **Temporal persistence**: `MeanPull` controls how quickly wind returns to prevailing patterns

# Constructor Example
    # Three turbines with westerly prevailing winds
    wind_triple = WindDirTriple(
        Init=[270.0, 275.0, 270.0],           # Slightly different targets per turbine
        CholSig=0.5 * I(3) + 0.2 * ones(3,3), # 0.5° individual + 0.2° shared noise
        MeanPull=0.15                         # 15% reversion per time step
    )

# Usage
This struct is used with the `Direction_RW_with_Mean` mode in wind field configurations:
    wind = Wind(
        input\\_dir="RW\\_with\\_Mean",
        dir=wind\\_triple
    )

See also: [`Direction_RW_with_Mean`](@ref), [`Wind`](@ref)
"""
struct WindDirTriple
    Init::Vector{Float64}      # Mean direction (vector or scalar)
    CholSig::Matrix{Float64}   # Cholesky factor of covariance matrix (nT x nT)
    MeanPull::Float64          # Scalar mean reversion factor
end

"""
    create_unified_buffers(wf::WindFarm, rotor_points=50) -> UnifiedBuffers

Create a unified buffer struct containing all arrays needed by interpolateOPs! and setUpTmpWFAndRun!.

# Arguments
- `wf::WindFarm`: Wind farm object to determine buffer sizes
- `rotor_points`: Number of rotor discretization points for FLORIS buffers (defaults to 50)

# Returns
- `UnifiedBuffers`: Struct containing all pre-allocated buffers including FLORIS computation buffers

# Note
For optimal performance, use the version that accepts a Floris object to automatically 
determine the correct rotor discretization size.
"""
function create_unified_buffers(wf::WindFarm, rotor_points=50)
    # For interpolateOPs!
    dist_buffer = zeros(wf.nOP)
    sorted_indices_buffer = zeros(Int, wf.nOP)
    
    # For setUpTmpWFAndRun!
    nT_with_grid = wf.nT + 1  # Original turbines + 1 grid point
    max_deps = wf.nT + 1  # Grid point depends on all original turbines
    
    M_buffer = zeros(nT_with_grid, 3)
    iTWFState_buffer = zeros(size(wf.States_WF, 2))
    tmp_Tpos_buffer = zeros(max_deps, 3)
    tmp_WF_buffer = zeros(max_deps, size(wf.States_WF, 2))
    tmp_Tst_buffer = zeros(max_deps, size(wf.States_T, 2))
    dists_buffer = zeros(max_deps)
    plot_WF_buffer = zeros(max_deps, size(wf.States_WF, 2))
    plot_OP_buffer = zeros(max_deps, 2)
    
    # Create FLORIS buffers with specified number of rotor points
    n_floris_points = max(rotor_points, 1)
    
    # Try to create FLORISBuffers if available, otherwise use nothing
    floris_buffers = try
    FLORISBuffers(n_floris_points)
    catch
        nothing
    end
    
    # Prepare a WindFarm buffer for grid-point computations (GP)
    GP = deepcopy(wf)
    original_nT = wf.nT
    GP.nT = original_nT + 1
    GP.dep = Vector{Vector{Int64}}(undef, GP.nT)
    for i in 1:original_nT
        GP.dep[i] = Int64[]
    end
    GP.dep[end] = collect(1:original_nT)
    if !isempty(wf.StartI)
        GP.StartI = hcat(wf.StartI, [wf.StartI[end] + 1])
    else
        GP.StartI = reshape([1], 1, 1)
    end
    GP.posBase = vcat(wf.posBase, zeros(1, 3))
    GP.posNac = vcat(wf.posNac, zeros(1, 3))
    GP.States_T = vcat(wf.States_T, zeros(1, size(wf.States_T, 2)))
    GP.D = vcat(wf.D, [0.0])
    GP.intOPs = [zeros(length(GP.dep[iT]), 4) for iT in 1:GP.nT]

    return UnifiedBuffers(
        dist_buffer,
        sorted_indices_buffer,
        M_buffer,
        iTWFState_buffer,
        tmp_Tpos_buffer,
        tmp_WF_buffer,
        tmp_Tst_buffer,
        dists_buffer,
        plot_WF_buffer,
        plot_OP_buffer,
        floris_buffers,
        GP
    )
end

# Method dispatch for Floris objects - defined later after Floris is loaded
# This will be defined in floridyn_cl.jl after all includes are processed

include("visualisation/structs_measurements.jl")
include("settings.jl")

# functions for calculating the wind field
include("windfield/windfield_direction.jl")
include("windfield/windfield_shear.jl")
include("windfield/windfield_turbulence.jl")
include("windfield/windfield_velocity.jl")

include("floris/discretization.jl")
include("floris/gaussian.jl")
include("floris/runfloris.jl")
include("floridyn_cl/floridyn_cl.jl")

include("correction/direction.jl")
include("correction/velocity.jl")
include("correction/turbulence.jl")

include("floridyn_cl/prepare_simulation.jl")
include("floridyn_cl/iterate.jl")

include("controller/controller.jl")
include("visualisation/calc_flowfield.jl")
include("visualisation/calc_power.jl")
include("visualisation/plot_flowfield.jl")
include("visualisation/plot_measurements.jl")
include("visualisation/create_video.jl")
include("visualisation/high_res_time.jl")
include("visualisation/pretty_print.jl")
include("visualisation/smart_plotting.jl")

"""
    run_floridyn(plt, set, wf, wind, sim, con, vis, 
                 floridyn, floris; msr=VelReduction) -> (WindFarm, DataFrame, Matrix)

Unified function that automatically handles both multi-threading and single-threading modes
for running FLORIDyn simulations with appropriate plotting callbacks.

# Arguments
- `plt`: PyPlot instance, usually provided by ControlPlots
- `set`: Settings object. See: [Settings](@ref)
- `wf`: WindFarm struct. These are work arrays, not persistent objects. See: [WindFarm](@ref)
- `wind`: Wind field input settings. See: [Wind](@ref)
- `sim`: Simulation settings. See: [Sim](@ref)
- `con`: Controller settings. See: [Con](@ref)
- `vis`: Visualization settings. See: [Vis](@ref)
- `floridyn`: FLORIDyn model struct. See: [FloriDyn](@ref)
- `floris`: Floris model struct. See: [Floris](@ref)
- `msr`: Measurement index for online flow field plotting (VelReduction, AddedTurbulence or EffWind). 
         Default VelReduction. See: [MSR](@ref)

# Returns
- Tuple (wf, md, mi): WindFarm, measurement data, and interaction matrix
"""
function run_floridyn(plt, set, wf, wind, sim, con, vis, floridyn, floris; msr=VelReduction)
    if Threads.nthreads() > 1 && nprocs() > 1
        # Multi-threading mode: use remote plotting callback
        # The rmt_plot_flow_field function should be defined via remote_plotting.jl
        try
            return runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris; rmt_plot_fn=Main.rmt_plot_flow_field, msr)
        catch e
            if isa(e, UndefVarError)
                error("rmt_plot_flow_field function not found in Main scope. Make sure to include remote_plotting.jl and call init_plotting() first.")
            else
                rethrow(e)
            end
        end
    else
        # Single-threading mode: no plotting callback
        return runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris; msr)
    end
end

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
        Pkg.add(["LaTeXStrings", "Timers", "TerminalPager", "DistributedNext", "ControlPlots"])
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
        wind, sim, con, floris, floridyn, ta, tp = setup(settings_file)
        set = Settings(wind, sim, con)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        wf = initSimulation(wf, sim)
        runFLORIDyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)
    end

end
end