# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function correctVel(::Velocity_None, set, T, Wind, SimTime, paramFLORIS, tmpM)
    # Get data
    U, Wind = getDataVel(set, Wind, T, SimTime, tmpM, paramFLORIS)

    # Correct Velocity
   T.States_WF[T[:StartI], 1] = U

    return T, Wind
end

function getDataVel(set, Wind, T, SimTime, tmpM, paramFLORIS)
    # Initialize U
    U = nothing

    # Determine which input mode to use
    if Wind.input_vel == "I_and_I"
        U, Wind.vel = getWindSpeedT(set.vel_mode, Wind.vel, collect(1:T[:nT]), SimTime,
                                   T.States_WF[T[:StartI], 2], paramFLORIS.p_p)

        if (SimTime - Wind.vel.StartTime) > Wind.vel.WSE.Offset
            # Ufree = Ueff / reduction
            U = U ./ tmpM[:, 1]
        end

    elseif Wind.input_vel == "ZOH_wErrorCov"
        U = getWindSpeedT(T[:States_WF][T[:StartI], 1], Wind.vel.ColSig)
    
    elseif Wind.input_vel == "RW_with_Mean"
        U = getWindSpeedT(T[:States_WF][T[:StartI], 1], Wind.vel)

    elseif Wind.input_vel == "EnKF_InterpTurbine"
        U = getWindSpeedT_EnKF(Wind.vel, collect(1:T[:nT]), SimTime)

    elseif Wind.input_vel == "EnKF_RW"
        U =T.States_WF[T[:StartI], 1]
        U = U .+ transpose(randn(length(U)) * Wind.vel.CholSig)

    elseif Wind.input_vel == "EnKF_ZOH"
        U =T.States_WF[T[:StartI], 1]

    elseif Wind.input_vel == "CLC_weighted_ZOH"
        U = T.C_Vel *T.States_WF[:, 1]

    else
        U = getWindSpeedT(set.vel_mode, Wind.vel, collect(1:T[:nT]), SimTime)
    end

    return U, Wind
end

