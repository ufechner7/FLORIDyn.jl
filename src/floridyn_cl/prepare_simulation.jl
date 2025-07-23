# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    prepareSimulation(set::Settings, wind::Wind, con::Con, floridyn::FloriDyn, floris::Floris, turbProp, sim::Sim)

Prepares the simulation environment for a wind farm analysis using the provided settings and parameters.

# Arguments
- `set::Settings`: Simulation settings containing configuration options.
- `wind::Wind`: Wind conditions or wind field data. See: [`Wind`](@ref) 
- `con::Con`: Controller parameters of the turbines.  See: [`Con`](@ref)
- `floridyn::FloriDyn`: Parameters specific to the FLORIDyn model. See: [`FloriDyn`](@ref)
- `floris::Floris`: Parameters specific to the FLORIS model. See: [`Floris`](@ref)
- `turbProp`: Properties of the turbines involved in the simulation.
- `sim::Sim`: Simulation-specific parameters or state. See: [`Sim`](@ref)

# Arguments that get modified
- `wind`: Updated with wind velocity, direction, turbulence intensity, and shear profile.
- `con`: Updated with yaw data.
- `sim`: Updated with the number of simulation steps.
- `floris`: May include additional parameters for the FLORIS model.

# Returns
- Returns the tuple `(T, wind, sim, con, floris)` where:
  - `T`: Dictionary containing turbine states and positions.
  - `wind`: Updated wind conditions.
  - `sim`: Updated simulation parameters.
  - `con`: Updated controller parameters.
  - `floris`: Parameters for the FLORIS model.
"""
function prepareSimulation(set::Settings, wind::Wind, con::Con, floridyn::FloriDyn, floris::Floris, turbProp, sim::Sim)
    loadDataWarnings = String[]

    # ========== WIND: Velocity ==========
    vel_file_dir = sim.path_to_data

    nT = size(turbProp.Pos, 1)
    input_vel = wind.input_vel

    if input_vel == "I_and_I"
        wind.vel.WSE = WSEParameters(nT, sim.path_to_data, sim.TimeStep)
        wind.vel.TimePrev = sim.start_time
        wind.vel.start_time = sim.start_time
    # elseif input_vel == "Interpolation"
    #     try
    #         wind.Vel = CSV.read("WindVel.csv", DataFrame)
    #     catch
    #         push!(loadDataWarnings, "WindVel.csv not found.")
    #     end
    # elseif input_vel == "InterpTurbine"
    #     try
    #         wind.Vel = CSV.read("WindVelTurbine.csv", DataFrame)
    #     catch
    #         push!(loadDataWarnings, "WindVelTurbine.csv not found, default created.")
    #     end
    elseif input_vel == "Constant"
        try 
            path = joinpath(vel_file_dir, "WindVelConstant.csv")
            if !isfile(path)
                @error "WindVelConstant.csv not found in $vel_file_dir"
            end
            df = CSV.read(path, DataFrame; header=false)
            wind.vel = df[1,1]
        catch
            push!(loadDataWarnings, "WindVelConstant.csv not found!")
        end
    end
    # elseif input_vel in ["ZOH_wErrorCov", "RW_with_Mean", "Interpolation_wErrorCov", "InterpTurbine_wErrorCov", "Constant_wErrorCov"]
    #     wind.Vel.Data = CSV.read("WindVelConstant.csv", DataFrame)
    #     VelCov = CSV.read("WindVelCovariance.csv", DataFrame)
    #     # Assumes readCovMatrix returns tuple
    #     _, wind.Vel.CholSig = readCovMatrix(VelCov, nT, "WindVel")
    #     if input_vel == "RW_with_Mean"
    #         wind.Vel.MeanPull = 1.0
    #     end
    # else
    #     error("Unknown wind velocity method: $input_vel")
    # end

    # Assuming the following global-like variables/structures:
    # wind.input_dir : String
    # wind.dir : Custom structure or Dict
    # sim.path_to_data : String
    # sim.start_time : Float64 (or other numeric)
    # sim.end_time : Float64
    # turbProp.Pos : Matrix or Vector
    # loadDataWarnings : Vector{String}
    # Define your own readCovMatrix function before using them

    nT = size(turbProp.Pos, 1)
    data_path = sim.path_to_data

    if wind.input_dir == "Interpolation"
        try
            path = joinpath(data_path, "WindDir.csv")
            wind.dir = readdlm(path, ',', Float64)
        catch
            push!(loadDataWarnings, "WindDir.csv not found.")
        end
    elseif wind.input_dir == "InterpTurbine"
        try
            path = joinpath(data_path, "WindDirTurbine.csv")
            wind.dir = readdlm(path, ',', Float64)
        catch
            push!(loadDataWarnings, "WindDirTurbine.csv not found")
        end
    elseif wind.input_dir == "Constant"
        try
            path = joinpath(data_path, "WindDirConstant.csv")
            wind.dir = readdlm(path, ',', Float64)
        catch
            push!(loadDataWarnings, "WindDirConstant.csv not found.")
        end
    elseif wind.input_dir == "Interpolation_wErrorCov"
        wind.dir = Dict()
        wind.dir[:Data] = readdlm("WindDir.csv", ',', Float64)
        DirCov = readdlm("WindDirCovariance.csv", ',', Float64)
        _, wind.dir[:CholSig] = readCovMatrix(DirCov, nT, "WindDir")
    elseif wind.input_dir == "InterpTurbine_wErrorCov"
        wind.dir = Dict()
        wind.dir[:Data] = readdlm("WindDirTurbine.csv", ',', Float64)
        DirCov = readdlm("WindDirCovariance.csv", ',', Float64)
        _, wind.dir[:CholSig] = readCovMatrix(DirCov, nT, "WindDir")
    elseif wind.input_dir == "Constant_wErrorCov"
        wind.dir = Dict()
        wind.dir[:Data] = readdlm("WindDirConstant.csv", ',', Float64)
        DirCov = readdlm("WindDirCovariance.csv", ',', Float64)
        _, wind.dir[:CholSig] = readCovMatrix(DirCov, nT, "WindDir")
    elseif wind.input_dir == "RW_with_Mean"
        wind.dir = Dict()
        wind.dir[:Init] = readdlm("WindDirConstant.csv", ',', Float64)
        DirCov = readdlm("WindDirCovariance.csv", ',', Float64)
        _, wind.dir[:CholSig] = readCovMatrix(DirCov, nT, "WindDir")
        wind.dir[:MeanPull] = 1
    else
        error("Method for wind direction $(wind.input_dir) unknown.")
    end

    # ============= TI =============        
    if wind.input_ti == "Interpolation"
        try
            wind.ti = CSV.read("WindTI.csv", DataFrame)
        catch e
            push!(loadDataWarnings, "WindTI.csv not found.")
        end
    elseif wind.input_ti == "InterpTurbine"
        try
            path = joinpath(data_path, "WindTITurbine.csv")
            wind.ti = CSV.read(path, DataFrame)
        catch e
            push!(loadDataWarnings, "WindTITurbine.csv not found.")
        end
    elseif wind.input_ti == "Constant"
        try
            path = joinpath(data_path, "WindTIConstant.csv")
            df = CSV.read(path, DataFrame; header=false)
            wind.ti = df[1,1]
        catch e
            push!(loadDataWarnings, "WindTIConstant.csv not found.")
        end
    else
        error("Method for turbulence intensity $(wind.input_ti) unknown.")
    end

    # ============= WindShear =============
    if wind.input_shear == "PowerLaw"
        path = joinpath(data_path, "WindShearPowerLaw.csv")
        alpha = CSV.read(path, DataFrame; header=false)[1,1] # Assuming alpha is in the first row, first column
        z0 = 1.0 # Default roughness length
        wind.shear = WindShear(alpha, z0)
    elseif wind.input_shear == "Interpolation"
        path = joinpath(data_path, "WindShearProfile.csv")
        wind.shear = CSV.read(path, DataFrame)

    elseif wind.input_shear == "LogLaw"
        path = joinpath(data_path, "WindShearLogLaw.csv")
        wind.shear = CSV.read(path, DataFrame)
    else
        error("Method for wind shear $(wind.input_shear) unknown.")
    end

    # ========== Wind Farm Setup ==========
    T = WindFarm()
    T.posBase = turbProp.Pos
    T.nT = size(turbProp.Pos, 1)
        
        t_data = getTurbineData(turbProp.Type)
    T.posNac = t_data.NacPos
    T.D = t_data.D

    states = States()
    if floridyn.twf_model == "heterogeneous"
        push!(states.WF_names, "OP_ori")
        states.WF = length(states.WF_names)
    elseif floridyn.twf_model != "homogeneous"
        error("Unknown TWF model $(floridyn.twf_model). Use 'homogeneous' or 'heterogeneous'")
    end
    
    # OP State and turbine initialization
    n_op = floridyn.n_op
    T.States_OP = zeros(n_op *T.nT, states.OP)
    T.Names_OP = states.OP_names
    T.States_T  = zeros(n_op *T.nT, states.Turbine)
    T.Names_T   = states.T_names
    T.States_WF = zeros(n_op *T.nT, states.WF)
    T.Names_WF  = states.WF_names
    T.StartI    = collect(1:n_op:(n_op *T.nT))'
    T.nOP       = n_op
    T.red_arr   = ones(T.nT, T.nT)

    # # deltaUW fallback
    # if !haskey(floridyn, :deltaUW)
    #     floridyn.deltaUW = floridyn.deltaDW
    # end

    # ========== Control Setup ==========
    yaw_method = con.yaw
    if yaw_method == "Constant"
        try
            con.YawData = CSV.read("Control_YawConstant.csv", DataFrame)
        catch
            push!(loadDataWarnings, "Control_YawConstant.csv not found.")
        end
    elseif yaw_method == "InterpTurbine"
        try
            con.YawData = CSV.read("Control_YawInterpolation.csv", DataFrame)
        catch
            push!(loadDataWarnings, "Control_YawInterpolation.csv not found.")
        end
    elseif yaw_method == "SOWFA"
        nacelleYaw = importSOWFAFile(joinpath(vel_file_dir, "SOWFA_nacelleYaw.csv"))
        con.yaw_data = condenseSOWFAYaw([nacelleYaw[1:T.nT:end, 2] reshape(nacelleYaw[:,3],T.nT, :)'])
    else
        error("Unknown yaw method: $yaw_method")
    end

    # if !haskey(con, :tanhYaw)
    #     con.tanhYaw = false
    # end

    # # ========== Init State ===========
   T.States_OP, T.States_T, T.States_WF = InitStates(set, T, wind, turbProp.Init_States, floris, sim)

    # # ========== Simulation Setup ==========
    sim.n_sim_steps = length(sim.start_time:sim.time_step:sim.end_time)
    floris.rotor_points = sim.rotor_points

    # # ========== Visualization ==========
    # if Vis.FlowField.Plot.Online
    #     Vis.Film.MovFileEffU = joinpath(sim.path_to_data, "Results", "EffectiveWindSpeed.avi")
    #     Vis.Film.FrmFileEffU = Vector{Any}(undef, sim.nSimSteps)
    #     Vis.Film.InProgress = true
    # else
    #     Vis.Film.InProgress = false
    # end

    # if Vis.FlowField.Error.Online
    #     dirSnap = readdir(Vis.FlowField.Error.ValidationPath)
    #     steps = Int[]
    #     for name in dirSnap
    #         step = tryparse(Int, name)
    #         push!(steps, something(step, -1))
    #     end
    #     Vis.FlowField.Error.Steps = steps
    # end

    # # ========= Warnings if data was generated ==========
    if !isempty(loadDataWarnings)
        for w in loadDataWarnings
            @warn w
        end
        error("Data not loaded properly. Please provide the required files.")
    end

    return T, wind, sim, con, floris
end
