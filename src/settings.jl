# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

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
- `dir::Union{Nothing, Matrix{Float64}}`: Optional wind direction matrix.
- `ti::Union{Nothing, Float64}`: Optional turbulence intensity value.
- `shear::Union{Nothing, WindShear}`: Optional wind shear profile.
"""
@with_kw mutable struct Wind
    input_vel::String
    input_dir::String
    input_ti::String
    input_shear::String
    correction::WindCorrection
    pertubation::WindPerturbation
    vel::Union{Nothing, Float64} = nothing
    dir::Union{Nothing, Matrix{Float64}} = nothing
    ti::Union{Nothing, Float64} = nothing
    shear::Union{Nothing, WindShear} = nothing
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
    setup(filename)

Initializes or configures the system using the provided `filename`. The `filename` should specify the path to a configuration or settings file required for setup.

# Arguments
- `filename::String`: Path to the `.yaml` file to be used for setup.

# Returns
- The tuple `(wind, sim, con)` where:
  - `wind`: An instance of the [`Wind`](@ref) struct containing wind-related parameters.
  - `sim`: An instance of the [`Sim`](@ref) struct containing simulation parameters.
  - `con`: An instance of the [`Con`](@ref) struct containing controller parameters.
  - `floris`: An instance of the [`Floris`](@ref) struct containing FLORIS model parameters.
  - `floridyn`: An instance of the [`FloriDyn`](@ref) struct containing FLORIDyn model parameters.
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
    wind, sim, con, floris, floridyn
end

"""
    Settings(wind::Wind, sim::Sim, con::Con)

Create and return a [`Settings`](@ref) object using the provided `wind` and `sim` parameters.

# Arguments
- `wind::Wind`: An instance of the [`Wind`](@ref) struct containing wind-related parameters.
- `sim::Sim`: An instance of the [`Sim`](@ref) struct containing the simulation parameters.
- `con::Con`: An instance of the [`Con`](@ref) struct containing the controller parameters.

# Returns
- A `Settings` struct configured with the given wind and simulation parameters.

# Notes
- The function uses the `str2type` helper to convert string representations of model types 
  into their corresponding Julia types.
- The `Settings` struct encapsulates the model settings for velocity, direction, 
  turbulence intensity, shear, and correction modes.
"""
function Settings(wind::Wind, sim::Sim, con::Con)
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
             iterate_mode, control_mode)
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
    importSOWFAFile(filename::String, dataLines = 2:typemax(Int))

Reads from a custom-formatted text file and extracts columns: Turbine, Times, and nacelle.

- `filename`: The path to the input file.
- `dataLines`: A single range (e.g., `2:Inf`) or vector of tuple ranges (e.g., `[(2, Inf)]`) for rows to import.

Returns:
- A Matrix{Float64} with the selected column data.
"""
function importSOWFAFile(filename::String, dataLines = 2:typemax(Int))

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
    if dataLines isa UnitRange && dataLines != 2:typemax(Int)
        df = df[dataLines, :]
    elseif dataLines isa Vector
        keep_rows = falses(nrow(df))
        for (start, stop) in dataLines
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
