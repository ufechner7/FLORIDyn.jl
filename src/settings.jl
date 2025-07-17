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

function setup(filename)
    data = YAML.load_file(filename)
    wind_data = data["wind"]
    wind = convertdict(Wind, wind_data)
    sim_data = data["sim"]
    sim = convertdict(Sim, sim_data)
    con_data = data["con"]
    con = convertdict(Con, con_data)
    wind, sim, con
end

function Settings(wind::Wind, sim)
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

