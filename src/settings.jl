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

function setup(filename)
    data = YAML.load_file(filename)
    wind_data = data["wind"]
    wind = convertdict(Wind, wind_data)
end

