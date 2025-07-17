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
    iter_sigma_dw::Int
    iter_sigma_cw::Int
    iter_sigma_time::Int
end

@with_kw struct Dir
    iter_sigma_dw::Int
    iter_sigma_cw::Int
    iter_sigma_time::Int
end

@with_kw struct Dyn
    advection::Int
    advection_mod::String
    op_iteration::String
    op_iter_weights::Vector{Float64}
    vel::Vel
    dir::Dir
end

@with_kw struct Sim
    floris::String
    start_time::Int
    end_time::Int
    time_step::Int
    rotor_discret::String
    rotor_points::Int
    dyn::Dyn
    init::String
    save_init_state::Bool
    save_final_state::Bool
end

function setup(filename)
    data = YAML.load_file(filename)
    wind_data = data["wind"]
    wind = convertdict(Wind, wind_data)
    sim_data = data["sim"]
    sim = convertdict(Sim, sim_data)
    wind, sim
end

