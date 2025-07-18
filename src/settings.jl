# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

@with_kw struct WindPerturbation
    vel::Float64
    vel_sigma::Float64
    dir::Float64
    dir_sigma::Float64
    ti::Float64
    ti_sigma::Float64
end

@with_kw struct WindCorrection
    vel::String
    dir::String
    ti::String
end

@with_kw struct Wind
    input_vel::String
    input_dir::String
    input_ti::String
    input_shear::String
    correction::WindCorrection
    pertubation::WindPerturbation
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

@with_kw struct Sim
    floris::String
    start_time::Int64
    end_time::Int64
    time_step::Int64
    rotor_discret::String
    rotor_points::Int64
    dyn::Dyn
    init::String
    save_init_state::Bool
    save_final_state::Bool
end

@with_kw struct Con
    yaw::String
end

@with_kw struct Floris
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
end

@with_kw struct FloryDyn
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
  - `wind`: An instance of the `Wind` struct containing wind-related parameters.
  - `sim`: An instance of the `Sim` struct containing simulation parameters.
  - `con`: An instance of the `Con` struct containing controller parameters.
  - `floris`: An instance of the `Floris` struct containing FLORIS model parameters.
  - `florydyn`: An instance of the `FloryDyn` struct containing FLORIDyn model parameters.
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
    florydyn_data = data["florydyn"]
    florydyn = convertdict(FloryDyn, florydyn_data)
    wind, sim, con, floris, florydyn
end

"""
    Settings(wind::Wind, sim::Sim)

Create and return a [`Settings`](@ref) object using the provided `wind` and `sim` parameters.

# Arguments
- `wind::Wind`: An instance of the `Wind` struct containing wind-related parameters.
- `sim`: An instance of the `Sim` struct containing the simulation parameters.

# Returns
- A `Settings` struct configured with the given wind and simulation parameters.

# Notes
- The function uses the `str2type` helper to convert string representations of model types 
  into their corresponding Julia types.
- The `Settings` struct encapsulates various model configurations for velocity, direction, 
  turbulence intensity, shear, and correction modes.
"""
function Settings(wind::Wind, sim::Sim)
    vel_mode = str2type("Velocity_" * wind.input_vel)
    dir_mode = str2type("Direction_" * wind.input_dir)
    turb_mode = str2type("TI_" * wind.input_ti)
    shear_mode = str2type("Shear_" * wind.input_shear)
    cor_dir_mode = str2type("Direction_" * wind.correction.dir)
    cor_vel_mode = str2type("Velocity_" * wind.correction.vel)
    cor_turb_mode = str2type("TI_" * wind.correction.ti)
    iterate_mode = str2type(sim.dyn.op_iteration)
    Settings(vel_mode, dir_mode, turb_mode, shear_mode, cor_dir_mode, cor_vel_mode, cor_turb_mode, iterate_mode)
end

