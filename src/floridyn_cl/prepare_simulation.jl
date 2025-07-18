# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function prepareSimulation(wind, con, paramFLORIDyn, paramFLORIS, turbProp, sim)
    loadDataWarnings = String[]

    # ========== WIND: Velocity ==========
    vel_file_dir = sim.path_to_data

    nT = size(turbProp.Pos, 1)
    input_vel = wind.input_vel

    if input_vel == "I_and_I"
        wind.vel.WSE = WSEParameters(nT, sim.PathToSim, sim.TimeStep)
        wind.vel.TimePrev = sim.StartTime
        wind.vel.StartTime = sim.StartTime
    # elseif input_vel == "Interpolation"
    #     try
    #         wind.Vel = CSV.read("WindVel.csv", DataFrame)
    #     catch
    #         generateDemoCSV(vel_file_dir, "WindVel.csv", 2, nothing, [sim.StartTime, 8], [sim.EndTime, 10])
    #         push!(loadDataWarnings, "WindVel.csv not found, default created. Units: [s, ms^-1]")
    #     end
    # elseif input_vel == "InterpTurbine"
    #     try
    #         wind.Vel = CSV.read("WindVelTurbine.csv", DataFrame)
    #     catch
    #         generateDemoCSV(vel_file_dir, "WindVelTurbine.csv", 3, size(turbProp.Pos, 1), [sim.StartTime, 8], [sim.EndTime, 10])
    #         push!(loadDataWarnings, "WindVelTurbine.csv not found, default created. Units: [s, ms^-1]")
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

    # # ========== Turbine Setup ==========
    # T = Dict()
    # T[:posBase] = turbProp.Pos
    # T[:nT] = size(turbProp.Pos, 1)
    
    # t_data = getTurbineData(turbProp.Type)
    # T[:posNac] = t_data.NacPos
    # T[:D] = t_data.D

    # states = States()
    # if paramFLORIDyn.twf_model == "heterogeneous"
    #     push!(states.WF_names, "OP_ori")
    #     states.WF = length(states.WF_names)
    # elseif paramFLORIDyn.twf_model != "homogeneous"
    #     error("Unknown TWF model $(paramFLORIDyn.twf_model). Use 'homogeneous' or 'heterogeneous'")
    # end
    
    # # OP State and turbine initialization
    # n_op = paramFLORIDyn.n_op
    # T[:States_OP] = zeros(n_op * T[:nT], states.OP)
    # T[:Names_OP] = states.OP_names
    # T[:States_T]  = zeros(n_op * T[:nT], states.Turbine)
    # T[:Names_T]   = states.T_names
    # T[:States_WF] = zeros(n_op * T[:nT], states.WF)
    # T[:Names_WF]  = states.WF_names
    # T[:StartI]    = 1:n_op:(n_op * T[:nT])
    # T[:nOP]       = n_op
    # T[:red_arr]   = ones(T[:nT], T[:nT])

    # # deltaUW fallback
    # if !haskey(paramFLORIDyn, :deltaUW)
    #     paramFLORIDyn.deltaUW = paramFLORIDyn.deltaDW
    # end

    # # ========== Control Setup ==========
    # yaw_method = con.Yaw
    # if yaw_method == "Constant"
    #     try
    #         con.YawData = CSV.read("Control_YawConstant.csv", DataFrame)
    #     catch
    #         generateDemoCSV(vel_file_dir, "Control_YawConstant.csv", 1, nothing, 270, nothing)
    #         push!(loadDataWarnings, "Control_YawConstant.csv not found, default created. Unit: [deg]")
    #     end
    # elseif yaw_method == "InterpTurbine"
    #     try
    #         con.YawData = CSV.read("Control_YawInterpolation.csv", DataFrame)
    #     catch
    #         generateDemoCSV(vel_file_dir, "Control_YawInterpolation.csv", 3, nT, [sim.StartTime, 270], [sim.EndTime, 250])
    #         push!(loadDataWarnings, "Control_YawInterpolation.csv not found, default created. Unit: [deg]")
    #     end
    # elseif yaw_method == "SOWFA"
    #     nacelleYaw = importSOWFAFile(joinpath(vel_file_dir, "SOWFA_nacelleYaw.csv"))
    #     con.YawData = condenseSOWFAYaw([nacelleYaw[1:T[:nT]:end, 2] reshape(nacelleYaw[:,3], T[:nT], :)'])
    # else
    #     error("Unknown yaw method: $yaw_method")
    # end

    # if !haskey(con, :tanhYaw)
    #     con.tanhYaw = false
    # end

    # # ========== Init State ===========
    # T[:States_OP], T[:States_T], T[:States_WF] = InitStates(T, wind, turbProp.Init_States, paramFLORIS, sim)

    # # ========== Simulation Setup ==========
    # sim.nSimSteps = length(sim.StartTime:sim.TimeStep:sim.EndTime)
    # paramFLORIS.RotorPoints = sim.RotorPoints

    # # ========== Visualization ==========
    # if Vis.FlowField.Plot.Online
    #     Vis.Film.MovFileEffU = joinpath(sim.PathToSim, "Results", "EffectiveWindSpeed.avi")
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
        error("Data not loaded properly. Default files generated. Please overwrite with appropriate data.")
    end
    T = nothing

    return T, wind, sim, con, paramFLORIS
end
