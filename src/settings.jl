# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    WindPerturbation

A mutable struct for configuring stochastic wind perturbations in wind farm simulations.

This struct controls whether perturbations are applied to different wind parameters and defines
the magnitude of these perturbations using standard deviations. Wind perturbations are used to 
model uncertainty in wind measurements or to perform sensitivity analysis.

# Fields
- `vel::Bool`: Enable/disable velocity perturbations. When `true`, random perturbations are applied to wind velocity.
- `vel_sigma::Float64`: Standard deviation for velocity perturbations [m/s]. Determines the magnitude of random variations added to the wind velocity.
- `dir::Bool`: Enable/disable direction perturbations. When `true`, random perturbations are applied to wind direction.
- `dir_sigma::Float64`: Standard deviation for direction perturbations [degrees]. Determines the magnitude of random variations added to the wind direction.
- `ti::Bool`: Enable/disable turbulence intensity perturbations. When `true`, random perturbations are applied to turbulence intensity.
- `ti_sigma::Float64`: Standard deviation for turbulence intensity perturbations [-]. Determines the magnitude of random variations added to the turbulence intensity.

# Notes
- Perturbations are typically applied as additive Gaussian noise with zero mean and the specified standard deviation
- The perturbation flags (`vel`, `dir`, `ti`) act as switches to enable or disable specific types of perturbations
- Setting a flag to `false` will disable perturbations for that parameter regardless of the sigma value
- Standard deviations should be positive values.
"""
@with_kw mutable struct WindPerturbation
    vel::Bool
    vel_sigma::Float64
    dir::Bool
    dir_sigma::Float64
    ti::Bool
    ti_sigma::Float64
end

@with_kw struct WindCorrection
    vel::String
    dir::String
    ti::String
end

"""
    Wind

A mutable struct representing wind settings.

# Fields
- `input_vel::String`: The type of wind velocity input, e.g., "Constant", "Interpolation".
- `input_dir::String`: The type of wind direction input, e.g., "Constant", "Interpolation".
- `input_ti::String`: The type of turbulence intensity input, e.g., "Constant", "Interpolation".
- `input_shear::String`: The type of wind shear input, e.g., "PowerLaw", "Interpolation".
- `correction::WindCorrection`: Settings for wind corrections.
- `perturbation::WindPerturbation`: Settings for wind perturbations.
- `vel::Union{Nothing, Float64}`: Optional wind velocity value.
- `dir::Union{Nothing, Matrix{Float64}, WindDirMatrix, WindDirType, WindDirTriple}`: Optional wind direction matrix or covariance data.
- `ti::Union{Nothing, Float64}`: Optional turbulence intensity value.
- `shear::Union{Nothing, WindShear}`: Optional wind shear profile.
"""
@with_kw mutable struct Wind
    input_vel::String
    input_dir::String
    input_ti::String
    input_shear::String
    correction::WindCorrection
    perturbation::WindPerturbation
    vel::Union{Nothing, Float64} = nothing
    dir::Union{Nothing, Matrix{Float64}, WindDirMatrix, WindDirType, WindDirTriple} = nothing
    ti::Union{Nothing, Float64, Matrix{Float64}} = nothing
    shear::Union{Nothing, WindShear, Matrix{Float64}} = nothing
end

@with_kw struct Vel
    iter_sigma_dw::Int64
    iter_sigma_cw::Int64
    iter_sigma_time::Int64
end

@with_kw struct Dir
    iter_sigma_dw::Int64
    iter_sigma_cw::Int64
    iter_sigma_time::Int64
end

@with_kw struct Dyn
    advection::Int64
    advection_mod::String
    op_iteration::String
    op_iter_weights::Vector{Float64}
    vel::Vel
    dir::Dir
end

"""
    Sim

A mutable struct representing the simulation settings.

# Fields
- `floris::String`: The name of the FLORIS model to be used.
- `start_time::Int64`: The start time of the simulation in seconds.
- `end_time::Int64`: The end time of the simulation in seconds.
- `time_step::Int64`: The time step for the simulation in seconds.
- `sim_step::Union{Nothing, Int64}`: Optional simulation step size.
- `n_sim_steps::Union{Nothing, Int64}`: Optional number of simulation steps.
- `rotor_discret::String`: The rotor discretization method, e.g., "Uniform", "Gaussian".
- `rotor_points::Int64`: The number of rotor points for discretization.
- `dyn::Dyn`: The dynamic settings for the simulation.
- `init::String`: The initialization method, e.g., "init", "load".
- `path_to_data::String`: The path to the directory where simulation data is stored.
- `save_init_state::Bool`: Whether to save the initial state of the simulation.
- `save_final_state::Bool`: Whether to save the final state of the simulation.
"""
@with_kw mutable struct Sim
    floris::String
    start_time::Int64
    end_time::Int64
    time_step::Int64
    sim_step::Union{Nothing, Int64} = nothing
    n_sim_steps::Union{Nothing, Int64} = nothing
    rotor_discret::String
    rotor_points::Int64
    dyn::Dyn
    init::String
    path_to_data::String
    save_init_state::Bool
    save_final_state::Bool
end

"""
    Con

A mutable struct for configuration settings.

# Fields
- `yaw::String`: The yaw control strategy, e.g., "Constant", "Interpolation".
- `yaw_data::Union{Nothing, Matrix{Float64}}`: Optional yaw data matrix.
- `tanh_yaw::Bool`: Whether to use hyperbolic tangent yaw control.
"""
@with_kw mutable struct Con
    yaw::String
    yaw_data::Union{Nothing, Matrix{Float64}} = nothing
    tanh_yaw::Bool = false
end

"""
    Floris

A mutable struct representing the settings for the FLORIDyn simulation. 

# Fields
- `alpha::Float64`: The alpha parameter for the FLORIS model.
- `beta::Float64`: The beta parameter for the FLORIS model.
- `k_a::Float64`: The k_a parameter for the FLORIS model.
- `k_b::Float64`: The k_b parameter for the FLORIS model.
- `k_fa::Float64`: The k_fa parameter for the FLORIS model.
- `k_fb::Float64`: The k_fb parameter for the FLORIS model.
- `k_fc::Float64`: The k_fc parameter for the FLORIS model.
- `k_fd::Float64`: The k_fd parameter for the FLORIS model.
- `eta::Int`: The eta parameter for the FLORIS model.
- `p_p::Float64`: The p_p parameter for the FLORIS model.
- `airDen::Float64`: The air density for the FLORIS model.
- `TIexp::Int`: The turbulence intensity exponent for the FLORIS model.
- `rotor_points::Union{Nothing, Int64}`: Optional number of rotor points.
"""
@with_kw mutable struct Floris
    alpha::Float64
    beta::Float64
    k_a::Float64
    k_b::Float64
    k_fa::Float64
    k_fb::Float64
    k_fc::Float64
    k_fd::Float64
    eta::Int
    p_p::Float64
    airDen::Float64
    TIexp::Int
    rotor_points::Union{Nothing, Int64} = nothing
end

"""
    FloriDyn

A structure representing the settings for the FLORIDyn simulation environment.

# Fields
- `n_op::Int`: The number of operating points.
- `deltaUW::Float64`: The delta U wind speed perturbation.
- `deltaDW::Float64`: The delta D wind direction perturbation.
- `deltaCW::Float64`: The delta C wind turbulence intensity perturbation.
- `dynStateChange::String`: The type of dynamic state change, e.g., "Constant", "Interpolation".
- `twf_model::String`: The type of TWF (Turbine Wake Flow) model used, e.g., "Gaussian", "FLORIDyn".
"""
@with_kw struct FloriDyn
    n_op::Int
    deltaUW::Float64
    deltaDW::Float64
    deltaCW::Float64
    dynStateChange::String
    twf_model::String
end

"""
    TurbineArray

A structure representing the configuration and properties of a wind turbine array.

# Fields
- `pos::Matrix{Float64}`: A matrix containing the positions of turbines. Each row represents
                          a turbine with columns for x, y, and z coordinates (in meters).
- `type::Vector{String}`: A vector of strings specifying the type/model of each turbine.
- `init_States::Matrix{Float64}`: A matrix containing the initial states of each turbine.
                                  Each row represents a turbine with columns for:
                                  - Column 1: `a` - axial induction factor
                                  - Column 2: `yaw` - initial yaw angle (in degrees)
                                  - Column 3: `ti` - turbulence intensity

# Example
```julia
# Create a simple 2-turbine array
pos = [0.0 0.0 0.0; 500.0 0.0 0.0]  # Two turbines 500m apart
type = ["NREL_5MW", "NREL_5MW"]
init_states = [0.33 0.0 0.1; 0.33 0.0 0.1]  # Both start with same initial conditions
turbines = TurbineArray(pos, type, init_states)
```
"""
struct TurbineArray
    pos::Matrix{Float64}
    type::Vector{String}
    init_States::Matrix{Float64}
end

"""
    Vis

A mutable struct for visualization settings in wind farm simulations.

This struct controls visualization options during simulation runtime, particularly for
live plotting and animation features. It provides both stored fields and computed 
properties for flexible path management.

# Fields
- `online::Bool`: Enable/disable online visualization during simulation. When `true`, 
                  live plots and animations are displayed during the simulation run.
                  When `false`, visualization is disabled for faster computation.
- `save::Bool`: Enable/disable saving of plots to disk. When `true`, plots are saved
                to disk in PNG format. When `false`, plots are only displayed (default: `false`).
- `print_filenames::Bool`: Enable/disable printing of saved filenames to console. 
                          When `true`, filenames of saved plots are printed (default: `false`).
- `video_folder::String`: Relative path for the video output directory. Used when `online=true`
                         and `save=true` (default: `"video"`).
- `output_folder::String`: Relative path for the output directory. Used when `online=false`
                          and `save=true` (default: `"out"`).
- `unique_output_folder::Bool`: If `true`, a unique timestamped folder is created for each simulation run.
                              If `false`, files are saved directly to the output folder (default: `true`).
- `flow_fields::Vector{FlowField}`: List of flow field visualizations to be created. Each [`FlowField`](@ref) 
                                   contains configuration for name, online display, video creation, and skip control.
                                   See [`FlowField`](@ref) for details (default: `FlowField[]`).
- `measurements::Vector{Measurement}`: List of measurement visualizations to be created. Each [`Measurement`](@ref) 
                                      contains configuration for name, separated plotting, and skip control.
                                      See [`Measurement`](@ref) for details (default: `Measurement[]`).
- `v_min::Float64`: Minimum velocity value for color scale in effective wind speed 
                   visualizations (msr=EffWind). Used to set consistent color scale limits 
                   across animation frames (default: `2.0`).
- `v_max::Float64`: Maximum velocity value for color scale in effective wind speed 
                   visualizations (msr=EffWind). Used to set consistent upper limits across 
                   animation frames (default: `10.0`).
- `rel_v_min::Float64`: Minimum relative velocity value for velocity reduction 
                       visualizations (msr=VelReduction). Controls the color scale for relative 
                       wind speed plots. Range: \\[0, 100\\] (default: `20.0`).
- `rel_v_max::Float64`: Maximum relative velocity value for velocity reduction 
                       visualizations (msr=VelReduction). Controls the upper limit for relative 
                       wind speed plots. Range: \\[0, 100\\] (default: `100.0`).
- `turb_max::Float64`: Maximum turbulence value for added turbulence visualizations (msr=AddedTurbulence).
                      Controls the upper limit for turbulence plots (default: `35.0`).
- `up_int::Int`: Update interval - controls how frequently visualization updates occur.
                Higher values result in less frequent updates for better performance (default: `1`).
- `unit_test::Bool`: Enable unit test mode for visualization functions. When `true`, 
                    plots are automatically closed after 1 second for testing purposes (default: `false`).

# Computed Properties
- `video_path::String`: Full absolute path to the video directory. Automatically creates the 
                       directory if it doesn't exist. Path is environment-dependent: 
                       `~/scratch/video_folder` on Delft Blue, `pwd()/video_folder` elsewhere.
- `output_path::String`: Full absolute path to the output directory. Automatically creates the 
                        directory if it doesn't exist. Path is environment-dependent:
                        `~/scratch/output_folder` on Delft Blue, `pwd()/output_folder` elsewhere.

# Constructors
```julia
# Default constructor with keyword arguments
vis = Vis(online=true, save=true, v_min=2.0, v_max=12.0)

# Constructor from YAML file
vis = Vis("path/to/config.yaml")  # Loads from data["vis"] section
```

# Examples
```julia
# Enable online visualization with plot saving and custom color scales
vis = Vis(online=true, save=true, v_min=2.0, v_max=12.0, rel_v_min=20.0, rel_v_max=100.0, up_int=5)

# Display only, no saving with default color scales
vis = Vis(online=true, save=false, v_min=2.0, rel_v_min=20.0)

# Specify flow field and measurement visualizations to create
flow_fields = [FlowField("flow_field_vel_reduction", true, true), FlowField("flow_field_added_turbulence", false, false)]
measurements = [Measurement("msr_vel_reduction", true), Measurement("msr_added_turbulence", false)]
vis = Vis(online=true, save=true, flow_fields=flow_fields, measurements=measurements)

# Disable online visualization for batch processing
vis = Vis(online=false, save=false)

# Load from YAML configuration (automatically parses flow_fields and measurements)
vis = Vis("data/vis_default.yaml")

# Access computed properties (creates directories automatically)
println("Saving plots to: ", vis.video_path)  # When online=true
println("Output directory: ", vis.output_path) # When online=false
```

# YAML Configuration Format
When loading from YAML files, the flow_fields and measurements are automatically parsed:
```yaml
vis:
  online: true
  save: true
  flow_fields:
    - name: "flow_field_vel_reduction"
      online: true
      create_video: true
      skip: false
    - "flow_field_added_turbulence"  # Simple format, defaults applied
  measurements:
    - name: "msr_vel_reduction" 
      separated: true
      skip: false
    - "msr_added_turbulence"  # Simple format, defaults applied
```

# File Saving Behavior
When `save=true`, plots are saved as PNG files with descriptive names:
- `ff_velocity_reduction.png` - for velocity reduction plots (msr=VelReduction)
- `ff_added_turbulence.png` - for turbulence intensity plots (msr=AddedTurbulence)  
- `ff_wind_speed.png` - for effective wind speed plots (msr=EffWind)
- Time-stamped versions: `ff_velocity_reduction_t0120s.png` when time parameter is provided

Save location depends on the `online` setting:
- `online=true`: Files saved to `video_path` (for animations)
- `online=false`: Files saved to `output_path` (for final results)

# Performance Notes
- Online visualization significantly slows down simulation performance
- Useful for debugging, monitoring simulation progress, or creating videos of the simulation
- When disabled, visualization functions are skipped to improve computational efficiency
- `up_int` can be used to reduce visualization frequency and improve simulation speed
- Directory creation is automatic but occurs only when computed properties are accessed
- Flow fields and measurements with `skip=true` are ignored completely for optimal performance

# Environment Adaptation
The struct automatically adapts to different computing environments:
- **Delft Blue supercomputer**: Uses `~/scratch/` directory for ample storage space
- **Local systems**: Uses current working directory (`pwd()`)
- Detection is automatic via the [`isdelftblue()`](@ref) function

# See Also
- [`FlowField`](@ref): Configuration struct for flow field visualizations
- [`Measurement`](@ref): Configuration struct for measurement visualizations
- [`parse_flow_fields`](@ref): Function to convert YAML flow field configurations
- [`parse_measurements`](@ref): Function to convert YAML measurement configurations
"""
@with_kw mutable struct Vis
    online::Bool = false              # Enable/disable online visualization during simulation; temporary variable
    no_plots::Int64 = 0               # Number of plots created, used for reporting
    no_videos::Int64 = 0              # Number of videos created, used for reporting
    show_plots::Bool = true           # master switch: if false, don't show plots, but results are still saved
    save::Bool = false                # save plots to video or output folder
    save_results::Bool = false        # save simulation results as .jld2 files
    print_filenames::Bool = false     # if true, print the names of the saved files
    log_debug::Bool = false           # if true, enable debug level logging output
    unique_folder::String = ""        # this will be set when starting the simulation
    video_folder::String = "video"    # relative video folder path
    output_folder::String = "out"     # relative output folder path
    unique_output_folder::Bool = true # if true, for each simulation run a new folder is created
    skip_flow_fields::Bool = false    # if true, completely skip creation of flow field visualizations (overrides individual entries)
    skip_measurements::Bool = false   # if true, completely skip creation of measurement visualizations (overrides individual entries)
    field_limits_min::Vector{Float64} = [0.0, 0.0, 0.0]          # [xmin, ymin, zmin] in meters
    field_limits_max::Vector{Float64} = [3000.0, 3000.0, 400.0]  # [xmax, ymax, zmax] in meters
    field_resolution::Float64 = 20.0                             # Resolution of the field in meters
    flow_fields::Vector{FlowField} = FlowField[]  # list of flow field visualizations to create
    measurements::Vector{Measurement} = Measurement[]  # list of measurement visualizations to create (parsed from YAML)
    v_min::Float64 = 2
    v_max::Float64 = 10
    rel_v_min::Float64 = 20
    rel_v_max::Float64 = 100
    turb_max::Float64 = 35
    up_int::Int = 1  # update interval
    unit_test::Bool = false  # enable unit test mode for visualization functions
end

# Helper to resolve a file under local or package data directories
function _resolve_data_path(filename::String)
    # If the provided path already exists (absolute or relative), use it
    if isfile(filename)
        return filename
    end
    pkg_root = joinpath(dirname(pathof(@__MODULE__)), "..")
    candidates = [
        joinpath(pwd(), filename),                 # relative to CWD
        joinpath(pwd(), "data", filename),        # under local data/
        joinpath(pkg_root, "data", filename),     # under package data/
    ]
    for p in candidates
        if isfile(p)
            return p
        end
    end
    return filename  # fallback; YAML.load_file will throw a helpful error
end

# Constructor for Vis struct from YAML file
function Vis(filename::String)
    filename = _resolve_data_path(filename)
    data = YAML.load_file(filename)
    vis_data = data["vis"]
    
    # Extract flow_fields and measurements separately to avoid recursion
    flow_fields_raw = get(vis_data, "flow_fields", [])
    measurements_raw = get(vis_data, "measurements", [])
    
    # Remove flow_fields and measurements from vis_data to avoid conflicts during convertdict
    vis_data_cleaned = copy(vis_data)
    delete!(vis_data_cleaned, "flow_fields")
    delete!(vis_data_cleaned, "measurements")
    
    # Create Vis struct using convertdict for other fields
    vis = convertdict(Vis, vis_data_cleaned)
    
    # Manually set the flow_fields and measurements fields
    vis.flow_fields = parse_flow_fields(flow_fields_raw)
    vis.measurements = parse_measurements(measurements_raw)
    
    return vis
end

# Add computed properties for video_path and output_path
"""
    Base.getproperty(vis::Vis, name::Symbol)

Provides computed properties `video_path` and `output_path` for the `Vis` struct, which resolve the full paths for video and output folders, respectively.

If the requested property is `:video_path` or `:output_path`, this method constructs the corresponding path based on the current environment (using either the user's home directory or the current working directory) and ensures that the directory exists by creating it if necessary.

For all other properties, it returns the field value as usual.

# Side Effects
- Accessing `video_path` or `output_path` will automatically create the corresponding directory if it does not exist.

# Arguments
- `vis::Vis`: Visualization settings struct.
- `name::Symbol`: The property name to access.

# Returns
- The value of the requested property, or the computed path for `video_path`/`output_path`.
"""
function Base.getproperty(vis::Vis, name::Symbol)
    if name === :video_path
        # Refactored from ternary operator to explicit if/else for clarity
        if isdelftblue()
            path = joinpath(homedir(), "scratch", vis.video_folder, vis.unique_folder)
        else
            path = joinpath(pwd(), vis.video_folder, vis.unique_folder)
        end
        return String(rstrip(mkpath(path), ['/','\\']))
    elseif name === :output_path
        # Refactored from ternary operator to explicit if/else for clarity
        if isdelftblue()
            path = joinpath(homedir(), "scratch", vis.output_folder, vis.unique_folder)
        else
            path = joinpath(pwd(), vis.output_folder, vis.unique_folder)
        end
        return String(rstrip(mkpath(path), ['/','\\']))
    else
        return getfield(vis, name)
    end
end

"""
    setup(filename) -> (wind, sim, con, floris, floridyn, ta)

Load wind farm configuration from a YAML file and parse all simulation components.

This function reads a comprehensive wind farm configuration file and extracts all necessary 
parameters for setting up a FLORIDyn simulation, including wind conditions, simulation 
settings, control strategies, FLORIS model parameters, FLORIDyn dynamics, and turbine 
array layout.

# Arguments
- `filename::String`: Path to the YAML configuration file containing wind farm setup data.
  The file should contain sections for: wind, sim, con, floris, floridyn, and turbines.

# Returns
A 6-tuple containing fully configured simulation components:
- `wind::Wind`: Wind conditions and input specifications (velocity, direction, turbulence, shear, corrections, perturbations)
- `sim::Sim`: Simulation parameters (time range, discretization, dynamics, initialization, data paths)  
- `con::Con`: Control settings (yaw strategies, control data)
- `floris::Floris`: FLORIS wake model parameters (alpha, beta, k coefficients, air density, etc.)
- `floridyn::FloriDyn`: FLORIDyn dynamic model settings (operating points, perturbations, state changes)
- `ta::TurbineArray`: Turbine array configuration (positions, types, initial states)

# YAML File Structure
The configuration file must contain these top-level sections:
```yaml
wind:          # Wind input specifications and corrections
sim:           # Simulation time, discretization, and dynamics  
con:           # Control strategies and yaw data
floris:        # FLORIS model coefficients and parameters
floridyn:      # FLORIDyn dynamic model settings
turbines:      # Array of turbine definitions with position, type, initial states
```

# Example
```julia
# Load complete wind farm configuration
wind, sim, con, floris, floridyn, ta = setup("data/2021_9T_Data.yaml")

# Access turbine positions  
println("Number of turbines: ", size(ta.pos, 1))
println("Simulation duration: ", sim.end_time - sim.start_time, " seconds")
```

# See Also
- [`Wind`](@ref): Wind conditions and input specifications
- [`Sim`](@ref): Simulation parameters and settings  
- [`Con`](@ref): Control configuration
- [`Floris`](@ref): FLORIS wake model parameters
- [`FloriDyn`](@ref): FLORIDyn dynamic model settings
- [`TurbineArray`](@ref): Turbine array layout and properties
"""
function setup(filename)
    data = YAML.load_file(filename)
    wind_data = data["wind"]
    wind = convertdict(Wind, wind_data)
    sim_data = data["sim"]
    sim = convertdict(Sim, sim_data)
    con_data = data["con"]
    con = convertdict(Con, con_data)
    floris_data = data["floris"]
    floris = convertdict(Floris, floris_data)
    floridyn_data = data["floridyn"]
    floridyn = convertdict(FloriDyn, floridyn_data)
    turbines = data["turbines"]

    # Extract Position: x, y, z
    pos = [Float64[t["x"], t["y"], t["z"]] for t in turbines]
    pos = reduce(vcat, [p' for p in pos])  # transpose and concatenate into matrix (9×3)

    # Extract Type
    type = [String(t["type"]) for t in turbines]

    # Extract Init States: a, yaw, ti
    init_states = [Float64[t["a"], t["yaw"], t["ti"]] for t in turbines]
    init_states = reduce(vcat, [s' for s in init_states])  # transpose and concatenate

    ta = TurbineArray(pos, type, init_states)
    wind, sim, con, floris, floridyn, ta
end

"""
    Settings(wind::Wind, sim::Sim, con::Con, parallel=false, threading=false)

Create and return a [`Settings`](@ref) object using the provided `wind` and `sim` parameters.

# Arguments
- `wind::Wind`: An instance of the [`Wind`](@ref) struct containing wind-related parameters.
- `sim::Sim`: An instance of the [`Sim`](@ref) struct containing the simulation parameters.
- `con::Con`: An instance of the [`Con`](@ref) struct containing the controller parameters.
- `parallel::Bool`:  Enable plotting in a separate process (default: `false`)
- `threading::Bool`: Enable threading for parallel computation within a single process (default: `false`)

# Returns
- A `Settings` struct configured with the given wind and simulation parameters.

# Notes
- The function uses the `str2type` helper to convert string representations of model types 
  into their corresponding Julia types.
- The `Settings` struct encapsulates the model settings for velocity, direction, 
  turbulence intensity, shear, and correction modes.
"""
function Settings(wind::Wind, sim::Sim, con::Con, parallel=false, threading=false)
    vel_mode = str2type("Velocity_" * wind.input_vel)
    dir_mode = str2type("Direction_" * wind.input_dir)
    turb_mode = str2type("TI_" * wind.input_ti)
    shear_mode = str2type("Shear_" * wind.input_shear)
    cor_dir_mode = str2type("Direction_" * wind.correction.dir)
    cor_vel_mode = str2type("Velocity_" * wind.correction.vel)
    cor_turb_mode = str2type("TI_" * wind.correction.ti)
    iterate_mode = str2type(sim.dyn.op_iteration)
    control_mode = str2type("Yaw_" * con.yaw)
    Settings(vel_mode, dir_mode, turb_mode, shear_mode, cor_dir_mode, cor_vel_mode, cor_turb_mode, 
             iterate_mode, control_mode, parallel, threading)
end

"""
    getTurbineData(names::Vector{String}) -> NamedTuple

Retrieve nacelle positions and rotor diameters for a given list of wind turbine types.

# Arguments
- `names::Vector{String}`: A vector of wind turbine type names. Supported types include:
  - `"DTU 10MW"`
  - `"DTU 5MW"`
  - `"Senvion 6.2M"`
  - `"V116"`
  - `"V117"`
  - `"V162"`
  - `"GE Haliade X"`

# Returns
- A `NamedTuple` with the following fields:
  - `NacPos::Matrix{Float64}`: An `N × 3` matrix where each row corresponds to the (x, y, z) coordinates of the nacelle position for each turbine.
  - `D::Vector{Float64}`: A vector of rotor diameters corresponding to each turbine.

# Raises
- `ArgumentError` if an unknown or misspelled turbine name is encountered.

"""
function getTurbineData(names::Vector{String})
    num = length(names)
    # Initialize NacPos as a num x 3 Matrix
    NacPos = zeros(Float64, num, 3)
    D = zeros(Float64, num)

    for i in 1:num
        name = names[i]
        if name == "DTU 10MW"
            NacPos[i, :] = [0.0, 0.0, 119.0]
            D[i] = 178.4
        elseif name == "DTU 5MW"
            NacPos[i, :] = [0.0, 0.0, 119.0]  # Placeholder
            D[i] = 178.4                      # Placeholder
        elseif name == "Senvion 6.2M"
            NacPos[i, :] = [0.0, 0.0, 152.0 - 29.0]
            D[i] = 126.0
        elseif name == "V116"
            NacPos[i, :] = [0.0, 0.0, 84.0]
            D[i] = 116.0
        elseif name == "V117"
            NacPos[i, :] = [0.0, 0.0, 84.0]
            D[i] = 117.0
        elseif name == "V162"
            NacPos[i, :] = [0.0, 0.0, 119.0]
            D[i] = 162.0
        elseif name == "GE Haliade X"
            NacPos[i, :] = [0.0, 0.0, 150.0]
            D[i] = 220.0
        else
            error("Turbine type '$name' not known or misspelled.")
        end
    end

    return (NacPos = NacPos, D = D)
end

"""
    importSOWFAFile(filename, data_lines = 2:typemax(Int))

Reads from a custom-formatted text file and extracts columns: Turbine, Times, and nacelle.

- `filename`: The path to the input file.
- `data_lines`: A single range (e.g., `2:Inf`) or vector of tuple ranges (e.g., `[(2, Inf)]`) for rows to import.

Returns:
- A Matrix{Float64} with the selected column data.
"""
function importSOWFAFile(filename, data_lines = 2:typemax(Int))
    pkg_path = joinpath(dirname(pathof(@__MODULE__)), "..")
    if ! isfile(filename)
        filename = joinpath(pkg_path, filename)
    end
    # Read full table first
    df = CSV.read(filename, DataFrame;
        delim=' ',
        ignorerepeated=true,
        missingstring="",
        header=[:Turbine, :Times, :Var3, :nacelle],
        types=Dict(:Turbine=>Float64, :Times=>Float64, :nacelle=>Float64, :Var3=>String),
        silencewarnings=true,
        ignoreemptyrows=true,
        skipto=2
    )

    # Filter rows if needed
    if data_lines isa UnitRange && data_lines != 2:typemax(Int)
        df = df[data_lines, :]
    elseif data_lines isa Vector
        keep_rows = falses(nrow(df))
        for (start, stop) in data_lines
            keep_rows[start:min(stop, nrow(df))] .= true
        end
        df = df[keep_rows, :]
    end

    # Select specific columns
    selected_df = df[:, [:Turbine, :Times, :nacelle]]

    # Convert to Matrix
    nacelleYaw = Matrix(selected_df)
    return nacelleYaw
end

function condenseSOWFAYaw(YawData::Array{wf,2}) where wf
    # Compute difference between adjacent rows (excluding the first and last rows)
    diff = sum(
        abs.(YawData[2:end-1, 2:end] .- YawData[1:end-2, 2:end]) .+
        abs.(YawData[2:end-1, 2:end] .- YawData[3:end, 2:end]),
        dims=2
    )
    # Find the rows where the difference is significant
    ind_important = vcat(1, findall(x -> x > 0, vec(diff)) .+ 1, size(YawData, 1))
    # Select the important rows
    return YawData[ind_important, :]
end

"""
    isdelftblue() -> Bool

Check if the current environment is the Delft Blue supercomputer.

This function determines whether the code is running on the Delft Blue supercomputer
by checking for the existence of the `~/scratch` directory, which is a characteristic
feature of the Delft Blue file system.

# Returns
- `Bool`: `true` if running on Delft Blue (i.e., `~/scratch` directory exists), 
          `false` otherwise.

# Usage
This function is used throughout FLORIDyn.jl to automatically adapt file paths 
to the computing environment. On Delft Blue, output files are stored in the 
scratch directory, which can be accessed remotely and offers lots of space.

# Example
```julia
if isdelftblue()
    output_path = joinpath(homedir(), "scratch", "out")
else
    output_path = joinpath(pwd(), "out")
end
```
"""
function isdelftblue()
    path = expanduser("~/scratch")
    ispath(path)
end

"""
        get_default_project() -> Tuple{String, String}

Read or create the default project selection and return `(settings_file, vis_file)` paths.

Behavior:
- Ensures `data/default.yaml` exists. If missing, it is created using the first project
    listed in `data/projects.yaml`.
- Reads `data/projects.yaml` and selects the project whose name matches `default.name`.
- Returns a tuple of full paths `(settings_file, vis_file)`, preferring files under the
    local `data/` folder and falling back to the package `data/` folder for reading.

Lookup strategy:
- Prefer files in the current working directory under `data/`.
- If not found, fall back to the package's `data/` directory for reading (not for writing).
"""
function get_default_project()
    # Resolve package root for read fallbacks
    pkg_root = joinpath(dirname(pathof(@__MODULE__)), "..")

    # Paths (prefer local workspace data dir)
    data_dir_local = joinpath(pwd(), "data")
    default_path_local = joinpath(data_dir_local, "default.yaml")
    projects_path_local = joinpath(data_dir_local, "projects.yaml")

    # Locate projects.yaml (read-only). Fall back to pkg data if not in local workspace.
    projects_path = projects_path_local
    if !isfile(projects_path)
        projects_path = joinpath(pkg_root, "data", "projects.yaml")
    end
    if !isfile(projects_path)
        error("projects.yaml not found in data directory (searched: $(projects_path_local))")
    end

    projects_data = YAML.load_file(projects_path)
    projects_list = get(projects_data, "projects", Any[])
    if isempty(projects_list)
        error("projects.yaml contains no projects")
    end

    # Helper to extract first project's name
    first_project = projects_list[1]["project"]
    first_name = String(first_project["name"])

    # Read or create default.yaml in the local workspace
    default_name = nothing
    default_msr = VelReduction  # Default fallback
    if isfile(default_path_local)
        try
            def_data = YAML.load_file(default_path_local)
            if haskey(def_data, "default") && haskey(def_data["default"], "name")
                default_name = String(def_data["default"]["name"])
            end
            if haskey(def_data, "default") && haskey(def_data["default"], "msr")
                msr_str = String(def_data["default"]["msr"])
                default_msr = toMSR(msr_str)
            end
        catch
            # If malformed, recreate from first project below
        end
    end
    if default_name === nothing
        # Create local data dir and write default.yaml with first project
        mkpath(data_dir_local)
        open(default_path_local, "w") do io
            write(io, "default:\n  name: $(first_name)\n  msr: $(string(default_msr))  # valid options: VelReduction, AddedTurbulence, EffWind\n")
        end
        default_name = first_name
    end

    # Find matching project by name
    chosen = nothing
    for entry in projects_list
        p = entry["project"]
        if String(p["name"]) == default_name
            chosen = p
            break
        end
    end
    if chosen === nothing
        # Fallback to first project and update default.yaml accordingly
        chosen = first_project
        open(default_path_local, "w") do io
            write(io, "default:\n  name: $(String(chosen["name"]))\n  msr: $(string(default_msr))  # valid options: VelReduction, AddedTurbulence, EffWind\n")
        end
    end

    # Build settings and vis file paths (prefer local workspace, fall back to package data)
    proj_name = String(chosen["name"])                  # e.g. "2021_9T_Data"
    vis_fname = String(chosen["vis"])                   # e.g. "vis_default.yaml"
    settings_fname = proj_name * ".yaml"                # e.g. "2021_9T_Data.yaml"

    # Candidate paths
    settings_local = joinpath(data_dir_local, settings_fname)
    settings_pkg   = joinpath(pkg_root, "data", settings_fname)
    vis_local      = joinpath(data_dir_local, vis_fname)
    vis_pkg        = joinpath(pkg_root, "data", vis_fname)

    # Resolve existing files
    settings_file = isfile(settings_local) ? settings_local : settings_pkg
    vis_file      = isfile(vis_local)      ? vis_local      : vis_pkg

    if !isfile(settings_file)
        error("Settings file not found: $(settings_fname) (searched in data/ and package data/)")
    end
    if !isfile(vis_file)
        error("Vis file not found: $(vis_fname) (searched in data/ and package data/)")
    end

    return (settings_file, vis_file)
end

"""
    list_projects() -> Vector{Tuple{String,String}}

Return a list of available projects as tuples `(name, vis)` using the same
projects.yaml discovery logic as `get_default_project()`.
"""
function list_projects()
    pkg_root = joinpath(dirname(pathof(@__MODULE__)), "..")
    projects_path_local = joinpath(pwd(), "data", "projects.yaml")
    projects_path = isfile(projects_path_local) ? projects_path_local : joinpath(pkg_root, "data", "projects.yaml")
    if !isfile(projects_path)
        error("projects.yaml not found in data directory (searched: $(projects_path_local))")
    end
    projects_data = YAML.load_file(projects_path)
    projects_list = get(projects_data, "projects", Any[])
    [(String(p["project"]["name"]), String(p["project"]["vis"])) for p in projects_list]
end

"""
    select_project() -> String

Interactive project selector. Tries a simple GUI (Gtk4.jl) if available, otherwise falls
back to a CLI prompt. Writes the chosen name to `data/default.yaml` and returns the
selected project name.

Usage:
    using FLORIDyn; select_project()
"""
function select_project()
    projs = list_projects()
    isempty(projs) && error("No projects found in projects.yaml")

    chosen_name::Union{Nothing,String} = nothing

    println("Available projects:")
    for (i,(n,_)) in enumerate(projs)
        println(rpad(string(i)*".",4), n)
    end
    print("Enter number [1-$(length(projs))]: ")
    flush(stdout)
    line = readline(stdin)
    idx = try parse(Int, strip(line)) catch; 0 end
    if idx < 1 || idx > length(projs)
        error("Invalid selection: $(line)")
    end
    chosen_name = projs[idx][1]

    # Write to local default.yaml
    data_dir_local = joinpath(pwd(), "data")
    mkpath(data_dir_local)
    default_path_local = joinpath(data_dir_local, "default.yaml")
    
    # Read existing MSR if available
    existing_msr = VelReduction  # Default fallback
    if isfile(default_path_local)
        try
            def_data = YAML.load_file(default_path_local)
            if haskey(def_data, "default") && haskey(def_data["default"], "msr")
                msr_str = String(def_data["default"]["msr"])
                existing_msr = toMSR(msr_str)
            end
        catch
            # If malformed, will use fallback
        end
    end
    
    open(default_path_local, "w") do io
        write(io, "default:\n  name: $(chosen_name)\n  msr: $(string(existing_msr))  # valid options: VelReduction, AddedTurbulence, EffWind\n")
    end
    println("Selected project saved to data/default.yaml: ", chosen_name)
    return chosen_name
end

"""
    get_default_msr() -> MSR

Read the default measurement type (MSR) from `data/default.yaml`.
If the file doesn't exist or doesn't have an `msr` field, returns `VelReduction`.

# Returns
- `MSR`: The default measurement type
"""
function get_default_msr()
    data_dir_local = joinpath(pwd(), "data")
    default_path_local = joinpath(data_dir_local, "default.yaml")
    
    if !isfile(default_path_local)
        return VelReduction  # Default fallback
    end
    
    try
        def_data = YAML.load_file(default_path_local)
        if haskey(def_data, "default") && haskey(def_data["default"], "msr")
            msr_str = String(def_data["default"]["msr"])
            return toMSR(msr_str)
        end
    catch
        # If malformed or missing, return default
    end
    
    return VelReduction  # Default fallback
end

"""
    set_default_msr(msr::MSR)

Set the default measurement type (MSR) in `data/default.yaml`.
Creates the file if it doesn't exist, preserving the existing project name.

# Arguments
- `msr::MSR`: The measurement type to set as default
"""
function set_default_msr(msr::MSR)
    data_dir_local = joinpath(pwd(), "data")
    default_path_local = joinpath(data_dir_local, "default.yaml")
    mkpath(data_dir_local)
    
    # Read existing data or use defaults
    default_name = "2021_54T_NordseeOne"  # fallback
    if isfile(default_path_local)
        try
            def_data = YAML.load_file(default_path_local)
            if haskey(def_data, "default") && haskey(def_data["default"], "name")
                default_name = String(def_data["default"]["name"])
            end
        catch
            # If malformed, will use fallback name
        end
    end
    
    # Write updated file
    open(default_path_local, "w") do io
        write(io, "default:\n  name: $(default_name)\n  msr: $(string(msr))  # valid options: VelReduction, AddedTurbulence, EffWind\n")
    end
    
    println("Default MSR saved to data/default.yaml: ", string(msr))
end

"""
    select_measurement() -> MSR

Interactive menu for selecting the default measurement type (MSR).
Displays available measurement types, prompts for user selection, and saves
the choice to `data/default.yaml`.

# Returns
- `MSR`: The selected measurement type

# Interactive Menu
The function displays:
1. VelReduction - Velocity reduction measurement
2. AddedTurbulence - Added turbulence measurement  
3. EffWind - Effective wind speed measurement

# Examples
```julia
# Run interactive selection
msr = select_measurement()
```

# Notes
- The selection is immediately saved to `data/default.yaml`
- The existing project name in `default.yaml` is preserved
- Invalid selections will throw an error and prompt retry
"""
function select_measurement()
    # Define measurement options with descriptions
    measurements = [
        (VelReduction, "VelReduction", "Velocity reduction measurement"),
        (AddedTurbulence, "AddedTurbulence", "Added turbulence measurement"),
        (EffWind, "EffWind", "Effective wind speed measurement")
    ]
    
    println()
    println("Available measurement types:")
    println("=" ^ 50)
    
    for (i, (msr, name, description)) in enumerate(measurements)
        println("$i. $name - $description")
    end
    
    println()
    print("Select measurement type (1-$(length(measurements))): ")
    line = readline()
    
    # Parse and validate selection
    idx = try parse(Int, strip(line)) catch; 0 end
    if idx < 1 || idx > length(measurements)
        error("Invalid selection: $(line). Please choose a number between 1 and $(length(measurements))")
    end
    
    selected_msr = measurements[idx][1]
    selected_name = measurements[idx][2]
    
    # Save the selection
    set_default_msr(selected_msr)
    
    println("Selected measurement type: ", selected_name)
    return selected_msr
end