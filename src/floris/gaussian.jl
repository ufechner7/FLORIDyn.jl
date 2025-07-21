# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    CalcCt(a, _)

Calculate the thrust coefficient (Ct) for a wind turbine based on the axial induction factor `a`.

# Arguments
- `a::Number`: Axial induction factor, typically between 0 and 0.5.
- _: unused parameter

# Returns
- `Ct::Number`: The calculated thrust coefficient.
"""
function CalcCt(a, _)
    Ct = 4 .* a .* (1 .- a)
    return Ct
end

mutable struct States
    T_names::Vector{String}
    Turbine::Int
    OP_names::Vector{String}
    OP::Int
    WF_names::Vector{String}
    WF::Int
end

function States()
    # Turbine states
    T_names = ["a", "yaw", "TI"]
    Turbine = length(T_names)

    # Observation point states
    OP_names = ["x0", "y0", "z0", "x1", "y1", "z1"]
    OP = length(OP_names)

    # Wind field states
    WF_names = ["wind_vel", "wind_dir", "TI0"]
    WF = length(WF_names)

    return States(T_names, Turbine, OP_names, OP, WF_names, WF)
end

"""
    Centerline(States_OP, States_T, States_WF, paramFLORIS, D)

Compute the centerline wake properties for a wind farm simulation.

# Arguments
- `States_OP`: Operational states of the turbines (e.g., yaw, pitch, etc.).
- `States_T`: Turbine-specific states (e.g., rotor speed, torque, etc.).
- `States_WF`: Wind farm-level states (e.g., wind direction, wind speed, etc.).
- `paramFLORIS`: Parameters for the FLORIS wake model.
- `D`: Rotor diameter or characteristic length scale.

# Returns
- The computed centerline wake properties `delta`, which includes the deflection in the y and z directions.

# Notes
This function is part of the Gaussian wake model implementation for wind farm simulations using the FLORIDyn.jl package.
"""
function Centerline(States_OP, States_T, States_WF, paramFLORIS, D)
    # Parameters
    k_a   = paramFLORIS.k_a
    k_b   = paramFLORIS.k_b
    alpha = paramFLORIS.alpha
    beta  = paramFLORIS.beta

    # States
    C_T   = CalcCt(States_T[:,1], States_T[:,2])
    yaw   = .-deg2rad.(States_T[:,2])
    I     = sqrt.(States_T[:,3].^2 .+ States_WF[:,3].^2)
    OPdw  = States_OP[:,4]

    # Calc x_0 (Core length)
    x_0 = (cos.(yaw) .* (1 .+ sqrt.(1 .- C_T)) ./ 
         (sqrt(2) .* (alpha .* I .+ beta .* (1 .- sqrt.(1 .- C_T))))) .* D

    # Calc k_z and k_y based on I
    k_y = k_a .* I .+ k_b
    k_z = k_y

    # Get field width y
    zs = zeros(size(OPdw))
    sig_y = max.(OPdw .- x_0, zs) .* k_y .+
        min.(OPdw ./ x_0, zs .+ 1) .* cos.(yaw) .* D ./ sqrt(8)

    # Get field width z
    sig_z = max.(OPdw .- x_0, zs) .* k_z .+
        min.(OPdw ./ x_0, zs .+ 1) .* D ./ sqrt(8)

    # Calc Theta
    Theta = 0.3 .* yaw ./ cos.(yaw) .* (1 .- sqrt.(1 .- C_T .* cos.(yaw)))

    # Calc Delta/Deflection
    delta_nfw = Theta .* min.(OPdw, x_0)

    delta_fw_1 = Theta ./ 14.7 .* sqrt.(cos.(yaw) ./ (k_y .* k_z .* C_T)) .* (2.9 .+ 1.3 .* sqrt.(1 .- C_T) .- C_T)
    delta_fw_2 = log.(Complex.(
        (1.6 .+ sqrt.(C_T)) .* 
        (1.6 .* sqrt.((8 .* sig_y .* sig_z) ./ (D^2 .* cos.(yaw))) .- sqrt.(C_T)) ./
        ((1.6 .- sqrt.(C_T)) .*
        (1.6 .* sqrt.((8 .* sig_y .* sig_z) ./ (D^2 .* cos.(yaw))) .+ sqrt.(C_T)))
    ))
    # println("delta_fw_2: ", delta_fw_2)
    # Use signbit and broadcasting for the sign/(2) + 0.5 logic
    factor = (sign.(OPdw .- x_0) ./ 2 .+ 0.5)

    deltaY = delta_nfw .+ factor .* delta_fw_1 .* delta_fw_2 .* D
    
    # Deflection in y and z direction
    delta = hcat(deltaY, zeros(size(deltaY)))
    
    return delta
end


function InitStates(set::Settings, T, Wind, InitTurb, paramFLORIS, Sim)
    # Unpack state arrays and parameters
    States_OP   = copy(T[:States_OP])
    States_T    = copy(T[:States_T])
    States_WF   = copy(T[:States_WF])
    nT          = T[:nT]
    nOP         = T[:nOP]
    deltaT      = Sim.time_step
    startTime   = Sim.start_time

    for iT = 1:nT
        # Retrieve wind field data
        if Wind.input_vel == "I_and_I"
            U = getWindSpeedT(set.vel_mode, Wind.vel, iT, startTime)
        elseif Wind.input_vel in ["ZOH_wErrorCov", "RW_with_Mean"]
            U = Wind.vel.Init
        else
            U = getWindSpeedT(set.vel_mode, Wind.vel, iT, startTime)
        end

        if Wind.input_dir == "RW_with_Mean"
            phiS = Wind.dir.Init
        else
            phiS = getWindDirT(set.dir_mode, Wind.dir, iT, startTime)
        end

        TI = getWindTiT(set.turb_mode, Wind.ti, iT, startTime)

        rangeOPs = ((iT-1)*nOP+1):(iT*nOP)

        # Initialize the States of the OPs and turbines
        States_WF[rangeOPs, 1] .= U
        States_WF[rangeOPs, 2] .= phiS
        States_WF[rangeOPs, 3] .= TI

        # Add orientation if used
        if length(T[:Names_WF]) == 4
            States_WF[rangeOPs, 4] .= phiS
        end

        # Downwind distance (wake coord)
        States_OP[rangeOPs, 4] .= (collect(0:(nOP-1)) .* deltaT .* U)

        # Init turbine states
        States_T[rangeOPs, :] = ones(nOP, 1) * InitTurb[iT, :]'

        # Crosswind position
        States_OP[rangeOPs, 5:6] = Centerline(States_OP[rangeOPs, :], States_T[rangeOPs, :],
                                              States_WF[rangeOPs, :], paramFLORIS, T[:D][iT])

        # Convert wind dir in fitting radians
        phiW = angSOWFA2world.(States_WF[rangeOPs, 2])

        # World coordinate position x0 and y0 including tower base and nacelle pos
        States_OP[rangeOPs, 1] .= cos.(phiW) .* States_OP[rangeOPs, 4] .-
                                   sin.(phiW) .* States_OP[rangeOPs, 5] .+
                                   T[:posBase][iT, 1] .+ T[:posNac][iT, 1]
        States_OP[rangeOPs, 2] .= sin.(phiW) .* States_OP[rangeOPs, 4] .+
                                   cos.(phiW) .* States_OP[rangeOPs, 5] .+
                                   T[:posBase][iT, 2] .+ T[:posNac][iT, 2]
        States_OP[rangeOPs, 3] .= States_OP[rangeOPs, 6] .+
                                   T[:posBase][iT, 3] .+ T[:posNac][iT, 3]
    end

    return States_OP, States_T, States_WF
end

function getVars(RPs, a, C_T, yaw, TI, TI0, param, D)
    # Unpack parameters
    k_a   = param.k_a
    k_b   = param.k_b
    alpha = param.alpha
    beta  = param.beta

    # States
    I = sqrt.(TI.^2 .+ TI0.^2)
    OPdw = RPs[:, 1]

    # Core length x_0
    x_0 = (cos.(yaw) .* (1 .+ sqrt.(1 .- C_T)) ./ 
          (sqrt(2) .* (alpha .* I .+ beta .* (1 .- sqrt.(1 .- C_T))))) .* D

    # Compute k_y and k_z
    k_y = k_a .* I .+ k_b
    k_z = k_y

    # Helper zero array
    zs = zeros(size(OPdw))

    # sig_y calculation (field width in y)
    sig_y = max.(OPdw .- x_0, zs) .* k_y .+
            min.(OPdw ./ x_0, zs .+ 1) .* cos.(yaw) .* D / sqrt(8)

    # sig_z calculation (field width in z)
    sig_z = max.(OPdw .- x_0, zs) .* k_z .+
            min.(OPdw ./ x_0, zs .+ 1) .* D / sqrt(8)

    # Theta
    Theta = 0.3 .* yaw ./ cos.(yaw) .* (1 .- sqrt.(1 .- C_T .* cos.(yaw)))

    # Deflection delta - near wake
    delta_nfw = Theta .* map((opdw, x0) -> min(opdw, x0), OPdw, x_0)

    # delta_fw parts
    delta_fw_1 = Theta ./ 14.7 .* sqrt.(cos.(yaw) ./ (k_y .* k_z .* C_T)) .* 
                 (2.9 .+ 1.3 .* sqrt.(1 .- C_T) .- C_T)

    # Intermediate term
    term = 1.6 .* sqrt.((8 .* sig_y .* sig_z) ./ (D.^2 .* cos.(yaw)))

    delta_fw_2 = log.(((1.6 .+ sqrt.(C_T)) .* (term .- sqrt.(C_T))) ./ 
                      ((1.6 .- sqrt.(C_T)) .* (term .+ sqrt.(C_T))))

    # Condition mask: OPdw > x_0 => 1.0, else 0.0
    mask = (OPdw .> x_0)
    blend = 0.5 .* sign.(OPdw .- x_0) .+ 0.5

    # Total delta in y
    deltaY = delta_nfw .+ blend .* delta_fw_1 .* delta_fw_2 .* D
    delta = hcat(deltaY, zeros(size(deltaY)))  # [delta_y, delta_z]

    # Potential core
    u_r_0 = (C_T .* cos.(yaw)) ./ 
            (2 .* (1 .- sqrt.(1 .- C_T .* cos.(yaw))) .* sqrt.(1 .- C_T))

    pc_y = D .* cos.(yaw) .* sqrt.(u_r_0) .* max.(1 .- OPdw ./ x_0, zs)
    pc_z = D .* sqrt.(u_r_0) .* max.(1 .- OPdw ./ x_0, zs)

    # For points exactly at the rotor plane
    rp = OPdw .== 0
    pc_y[rp] .= D .* cos.(yaw)
    pc_z[rp] .= D

    return sig_y, sig_z, C_T, x_0, delta, pc_y, pc_z
end

function runFLORIS(set::Settings, LocationT, States_WF, States_T, D, paramFLORIS, windshear)
    if D[end] > 0
        RPl, RPw = discretizeRotor(paramFLORIS.rotor_points)
    else
        RPl = [0.0 0.0 0.0]
        RPw = [1.0]
    end
    # Yaw rotation for last turbine
    tmp_yaw = deg2rad(States_T[end, 2])
    R = [cos(tmp_yaw)  sin(tmp_yaw)  0.0;
        -sin(tmp_yaw)  cos(tmp_yaw)  0.0;
         0.0           0.0           1.0]

    RPl = (R * (RPl .* D[end])')' .+ LocationT[end, :]'

    if length(D) == 1
        redShear = getWindShearT(set.shear_mode, windshear, RPl[:, 3] ./ LocationT[3])
        T_red_arr = RPw' * redShear
        T_aTI_arr, T_Ueff, T_weight = nothing, nothing, nothing
        return T_red_arr, T_aTI_arr, T_Ueff, T_weight
    end

    # Initialize outputs
    nT = length(D)
    T_red_arr = ones(nT)
    T_aTI_arr = zeros(nT - 1)
    T_weight = zeros(nT - 1)

    for iT in 1:(nT - 1)

        tmp_phi = size(States_WF,2) == 4 ? angSOWFA2world(States_WF[iT, 4]) :
                                           angSOWFA2world(States_WF[iT, 2])

        tmp_RPs = RPl .- LocationT[iT, :]'
        R_phi = [cos(tmp_phi)  sin(tmp_phi)  0.0;
                -sin(tmp_phi)  cos(tmp_phi)  0.0;
                 0.0           0.0           1.0]

        tmp_RPs = (R_phi * tmp_RPs')'

        if tmp_RPs[1, 1] <= 10
            continue
        end

        a = States_T[iT, 1]
        yaw_deg = States_T[iT, 2]
        yaw = -deg2rad(yaw_deg)
        TI = States_T[iT, 3]
        Ct = CalcCt(a, yaw_deg)
        TI0 = States_WF[iT, 3]

        sig_y, sig_z, C_T, x_0, delta, pc_y, pc_z = getVars(
            tmp_RPs, a, Ct, yaw, TI, TI0, paramFLORIS, D[iT]
        )

        cw_y = tmp_RPs[:, 2] .- delta[:, 1]
        cw_z = tmp_RPs[:, 3] .- delta[:, 2]
        phi_cw = atan.(cw_z, cw_y)
        r_cw = sqrt.(cw_y.^2 .+ cw_z.^2)

        core = (r_cw .< sqrt.((0.5 .* pc_y .* cos.(phi_cw)).^2 .+
                              (0.5 .* pc_z .* sin.(phi_cw)).^2)) .|
               (tmp_RPs[:, 1] .== 0.0)

        nw = tmp_RPs[:, 1] .< x_0

        tmp_RPs_r = zeros(size(RPw))
        tmp_RPs_r[core] .= 1 .- sqrt(1 .- C_T)

        fw = .!nw
        gaussAbs = zeros(size(RPw))
        gaussAbs[nw] .= 1 .- sqrt(1 .- C_T)
        gaussAbs[fw] .= 1 .- sqrt.(1 .- C_T .* cos(yaw) ./ (8 .* sig_y[fw] .* sig_z[fw] ./ D[iT]^2))

        gaussWght = ones(size(RPw))
        not_core = .!core
        if any(not_core)
            exp_y = @. exp(-0.5 * ((cw_y[not_core] - cos(phi_cw[not_core]) .* pc_y[not_core] * 0.5) ./ sig_y[not_core])^2)
            exp_z = @. exp(-0.5 * ((cw_z[not_core] - sin(phi_cw[not_core]) .* pc_z[not_core] * 0.5) ./ sig_z[not_core])^2)

            gaussWght[not_core] .= exp_y .* exp_z
            tmp_RPs_r[not_core] .= gaussAbs[not_core] .* gaussWght[not_core]
        end

        T_weight[iT] = sum(gaussWght)
        T_red_arr[iT] = 1 .- dot(RPw, tmp_RPs_r)

        # Added TI
        T_addedTI_tmp = paramFLORIS.k_fa * (
            a^paramFLORIS.k_fb *
            TI0^paramFLORIS.k_fc *
            (mean(tmp_RPs[:, 1]) / D[iT])^paramFLORIS.k_fd
        )

        TIexp = paramFLORIS.TIexp
        exp_y = @. exp(-0.5 * ((cw_y - cos(phi_cw) .* pc_y * 0.5) ./ (TIexp .* sig_y))^2)
        exp_z = @. exp(-0.5 * ((cw_z - sin(phi_cw) .* pc_z * 0.5) ./ (TIexp .* sig_z))^2)

        T_aTI_arr[iT] = T_addedTI_tmp * dot(RPw, exp_y .* exp_z)
    end

    redShear = getWindShearT(set.shear_mode, windshear, RPl[:, 3] ./ LocationT[end, 3])
    T_red_arr[end] = dot(RPw, redShear)

    T_red = prod(T_red_arr)
    T_Ueff = States_WF[end, 1] * T_red

    return T_red_arr, T_aTI_arr, T_Ueff, T_weight
end



