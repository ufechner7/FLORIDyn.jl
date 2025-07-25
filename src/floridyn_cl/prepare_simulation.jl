# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    readCovMatrix(cov_data, nT, name)

Read and process covariance matrix data for wind field error modeling.

# Arguments
- `cov_data`: Matrix containing covariance data (typically scalar or matrix)
- `nT::Int`: Number of turbines
- `name::String`: Name of the variable (for error messages)

# Returns
- Returns a tuple `(cov_matrix, chol_factor)` where:
  - `cov_matrix`: The full covariance matrix
  - `chol_factor`: Cholesky decomposition of the covariance matrix
"""
function readCovMatrix(cov_data, nT, name)
    if length(cov_data) == 1
        # Scalar variance - create diagonal covariance matrix
        var_val = cov_data[1]
        cov_matrix = var_val * I(nT)
        chol_factor = sqrt(var_val) * I(nT)
    else
        # Assume it's already a full covariance matrix
        cov_matrix = cov_data
        chol_factor = cholesky(cov_matrix).L
    end
    
    return cov_matrix, chol_factor
end

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
- Returns the tuple `(wf, wind, sim, con, floris)` where:
  - `wf`: Wind farm struct containing turbine states and positions. See: [`WindFarm`](@ref)
  - `wind`: Updated wind conditions.
  - `sim`: Updated simulation parameters.
  - `con`: Updated controller parameters.
  - `floris`: Parameters for the FLORIS model.
"""
function prepareSimulation(set::Settings, wind::Wind, con::Con, floridyn::FloriDyn, 
                           floris::Floris, turbProp, sim::Sim)
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
            # Generate demo CSV file for InterpTurbine mode
            generateDemoCSV(data_path, "WindDirTurbine.csv", 3, nT, [0.0, 270.0], [100.0, 270.0])
            wind.dir = readdlm(path, ',', Float64)
        end
    elseif wind.input_dir == "Constant"
        try
            path = joinpath(data_path, "WindDirConstant.csv")
            wind.dir = readdlm(path, ',', Float64)
        catch
            push!(loadDataWarnings, "WindDirConstant.csv not found.")
        end
    elseif wind.input_dir == "Interpolation_wErrorCov"
        data = readdlm(joinpath(data_path, "WindDir.csv"), ',', Float64)
        DirCov = readdlm(joinpath(data_path, "WindDirCovariance.csv"), ',', Float64)
        _, cholsig = readCovMatrix(DirCov, nT, "WindDir")
        wind.dir = WindDirMatrix(data, cholsig)
    elseif wind.input_dir == "InterpTurbine_wErrorCov"
        data = readdlm(joinpath(data_path, "WindDirTurbine.csv"), ',', Float64)
        DirCov = readdlm(joinpath(data_path, "WindDirCovariance.csv"), ',', Float64)
        _, cholsig = readCovMatrix(DirCov, nT, "WindDir")
        wind.dir = WindDirMatrix(data, cholsig)
    elseif wind.input_dir == "Constant_wErrorCov"
        data = readdlm(joinpath(data_path, "WindDirConstant.csv"), ',', Float64)
        DirCov = readdlm(joinpath(data_path, "WindDirCovariance.csv"), ',', Float64)
        _, cholsig = readCovMatrix(DirCov, nT, "WindDir")
        wind.dir = WindDirType(data[1], cholsig)
    elseif wind.input_dir == "RW_with_Mean"
        data = readdlm(joinpath(data_path, "WindDirConstant.csv"), ',', Float64)
        DirCov = readdlm(joinpath(data_path, "WindDirCovariance.csv"), ',', Float64)
        _, cholsig = readCovMatrix(DirCov, nT, "WindDir")
        init_values = fill(data[1], nT)  # Create vector of initial values for all turbines
        wind.dir = WindDirTriple(init_values, cholsig, 1.0)  # MeanPull = 1.0
    else
        error("Method for wind direction $(wind.input_dir) unknown.")
    end

    # ============= TI =============        
    if wind.input_ti == "Interpolation"
        path = joinpath(data_path, "WindTI.csv")
        try
            df = CSV.read(path, DataFrame)
            wind.ti = Matrix{Float64}(df)
        catch e
            push!(loadDataWarnings, "WindTI.csv not found.")
            generateDemoCSV(data_path, "WindTI.csv", 2, nT, [0.0, 0.0], [100.0, 100.0])
            df = CSV.read(path, DataFrame)
            wind.ti = Matrix{Float64}(df)
        end
    elseif wind.input_ti == "InterpTurbine"
        path = joinpath(data_path, "WindTITurbine.csv")
        try
            df = CSV.read(path, DataFrame)
            wind.ti = Matrix{Float64}(df)
        catch e
            push!(loadDataWarnings, "WindTITurbine.csv not found.")
            generateDemoCSV(data_path, "WindTITurbine.csv", 3, nT, [0.0, 0.06], [100.0, 0.06])
            df = CSV.read(path, DataFrame)
            wind.ti = Matrix{Float64}(df)
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
        df = CSV.read(path, DataFrame)
        wind.shear = Matrix{Float64}(df)

    elseif wind.input_shear == "LogLaw"
        path = joinpath(data_path, "WindShearLogLaw.csv")
        df = CSV.read(path, DataFrame)
        wind.shear = Matrix{Float64}(df)
    else
        error("Method for wind shear $(wind.input_shear) unknown.")
    end

    # ========== Wind Farm Setup ==========
    wf = WindFarm()
    wf.posBase = turbProp.Pos
    wf.nT = size(turbProp.Pos, 1)
        
        t_data = getTurbineData(turbProp.Type)
    wf.posNac = t_data.NacPos
    wf.D = t_data.D

    states = States()
    if floridyn.twf_model == "heterogeneous"
        push!(states.WF_names, "OP_ori")
        states.WF = length(states.WF_names)
    elseif floridyn.twf_model != "homogeneous"
        error("Unknown TWF model $(floridyn.twf_model). Use 'homogeneous' or 'heterogeneous'")
    end
    
    # OP State and turbine initialization
    n_op = floridyn.n_op
    wf.States_OP = zeros(n_op * wf.nT, states.OP)
    wf.Names_OP = states.OP_names
    wf.States_T  = zeros(n_op * wf.nT, states.Turbine)
    wf.Names_T   = states.T_names
    wf.States_WF = zeros(n_op * wf.nT, states.WF)
    wf.Names_WF  = states.WF_names
    wf.StartI    = collect(1:n_op:(n_op * wf.nT))'
    wf.nOP       = n_op
    wf.red_arr   = ones(wf.nT, wf.nT)

    # # deltaUW fallback
    # if !haskey(floridyn, :deltaUW)
    #     floridyn.deltaUW = floridyn.deltaDW
    # end

    # ========== Control Setup ==========
    yaw_method = con.yaw
    if yaw_method == "Constant"
        try
            df = CSV.read("Control_YawConstant.csv", DataFrame)
            con.yaw_data = Matrix{Float64}(df)
        catch
            push!(loadDataWarnings, "Control_YawConstant.csv not found.")
        end
    elseif yaw_method == "InterpTurbine"
        try
            df = CSV.read("Control_YawInterpolation.csv", DataFrame)
            con.yaw_data = Matrix{Float64}(df)
        catch
            push!(loadDataWarnings, "Control_YawInterpolation.csv not found.")
        end
    elseif yaw_method == "SOWFA"
        nacelleYaw = importSOWFAFile(joinpath(vel_file_dir, "SOWFA_nacelleYaw.csv"))
        con.yaw_data = condenseSOWFAYaw([nacelleYaw[1:wf.nT:end, 2] reshape(nacelleYaw[:,3],wf.nT, :)'])
    else
        error("Unknown yaw method: $yaw_method")
    end

    # if !haskey(con, :tanhYaw)
    #     con.tanhYaw = false
    # end

    # # ========== Init State ===========
   wf.States_OP, wf.States_T, wf.States_WF = init_states(set, wf, wind, turbProp.Init_States, floris, sim)

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
        # error("Data not loaded properly. Please provide the required files.")
    end

    return wf, wind, sim, con, floris
end

"""
    generateDemoCSV(path, name, type, nT, startV, endV)

Generate demo CSV files for wind field data based on the specified type.

# Arguments
- `path::String`: Directory path where the CSV file will be saved
- `name::String`: Name of the CSV file
- `type::Int`: Type of data generation (1=Constant, 2=Interpolation, 3=Turbine individual interpolation)
- `nT::Int`: Number of turbines
- `startV`: Starting value(s) - format depends on type
- `endV`: Ending value(s) - format depends on type

# Types
- Type 1 (Constant): `startV` should be a single value
- Type 2 (Interpolation): `startV` and `endV` should be [time, value] pairs
- Type 3 (Turbine individual interpolation): `startV` and `endV` should be [time, value] pairs

# Examples
```julia
# Constant value
generateDemoCSV("./data/", "WindDirConstant.csv", 1, 9, 270.0, nothing)

# Interpolation
generateDemoCSV("./data/", "WindDir.csv", 2, 9, [0.0, 250.0], [100.0, 280.0])

# Turbine individual interpolation  
generateDemoCSV("./data/", "WindDirTurbine.csv", 3, 9, [0.0, 250.0], [100.0, 280.0])
```
"""
function generateDemoCSV(path::String, name::String, type::Int, nT::Int, startV, endV)
    filepath = joinpath(path, name)
    
    if type == 1
        # Constant value
        # start: single value
        writedlm(filepath, startV, ',')
        
    elseif type == 2
        # Interpolation
        # start: [time, value]
        # end:   [time, value]
        data = vcat(startV', endV')
        @info "filepath: $filepath"
        @info "pwd(): $(pwd())"
        writedlm(filepath, data, ',')
        
    elseif type == 3
        # Turbine individual interpolation
        # start: [time, value]
        # end:   [time, value]
        A = [startV[2], endV[2]]
        T = [startV[1], endV[1]]
        data = hcat(T, repeat(A, 1, nT))
        writedlm(filepath, data, ',')
        
    else
        error("prepareSimulation -> generateDemoCSV: Type $type not defined.")
    end
    @info "Demo CSV file generated at: $filepath"
end
