# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function correctVel(::Velocity_None, set, wf, Wind, SimTime, paramFLORIS, tmpM)
    # Get data
    U, Wind = getDataVel(set, Wind, wf, SimTime, tmpM, paramFLORIS)

    # Correct Velocity
   wf.States_WF[wf.StartI, 1] = U

    return wf, Wind
end

function getDataVel(set, Wind, wf, SimTime, tmpM, paramFLORIS)
    # Initialize U
    U = nothing

    # Determine which input mode to use
    if Wind.input_vel == "I_and_I"
        U, Wind.vel = getWindSpeedT(set.vel_mode, Wind.vel, collect(1:wf.nT), SimTime,
                                   wf.States_WF[wf.StartI, 2], paramFLORIS.p_p)

        if (SimTime - Wind.vel.StartTime) > Wind.vel.WSE.Offset
            # Ufree = Ueff / reduction
            U = U ./ tmpM[:, 1]
        end

    elseif Wind.input_vel == "ZOH_wErrorCov"
        U = getWindSpeedT(wf.States_WF[wf.StartI, 1], Wind.vel.ColSig)
    
    elseif Wind.input_vel == "RW_with_Mean"
        U = getWindSpeedT(wf.States_WF[wf.StartI, 1], Wind.vel)

    elseif Wind.input_vel == "EnKF_InterpTurbine"
        U = getWindSpeedT_EnKF(Wind.vel, collect(1:wf.nT), SimTime)

    elseif Wind.input_vel == "EnKF_RW"
        U =wf.States_WF[wf.StartI, 1]
        U = U .+ transpose(randn(length(U)) * Wind.vel.CholSig)

    elseif Wind.input_vel == "EnKF_ZOH"
        U =wf.States_WF[wf.StartI, 1]

    elseif Wind.input_vel == "CLC_weighted_ZOH"
        U = wf.C_Vel *wf.States_WF[:, 1]

    else
        U = getWindSpeedT(set.vel_mode, Wind.vel, collect(1:wf.nT), SimTime)
    end

    return U, Wind
end

