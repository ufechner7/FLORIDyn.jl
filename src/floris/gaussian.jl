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


