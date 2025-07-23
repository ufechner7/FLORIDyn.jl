# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function correctTi(::TI_None, set, T, Wind, SimTime)
    # correctTi updates the turbulent intensity (TI) value inT.States_WF
    # at the rowT.StartI and column 3, using data from getDataTI.

    # Get new turbulent intensity value
    TI = getDataTI(set, Wind, T, SimTime)

    # Update the TI value in the States_WF matrix
    try
       T.States_WF[T[:StartI], 3] .= TI'
    catch e
        error("Error updatingT.States_WF: $(e.msg)")
    end

    return T
end

function getDataTI(set, Wind, T, SimTime)
    # GETDATATI retrieves the data for the ambient turbulence intensity

    if Wind.input_ti == "EnKF_InterpTurbine"
        TI = getWindTiT_EnKF(set.turb_mode, Wind.ti, collect(1:T.nT), SimTime)
    elseif Wind.input_ti == "EnKF_ZOH"
        TI =T.States_WF[T[:StartI], 3]
    elseif Wind.input_ti == "EnKF_RW"
        TI =T.States_WF[T[:StartI], 3]
        # 'phi' is not defined in your snippet. Replace with T.nT if that's correct.
        TI = TI .+ (randn(1, T.nT) * Wind.ti.CholSig)'
    elseif Wind.input_ti == "CLC_weighted_ZOH"
        TI = T.C_TI *T.States_WF[:, 3]
    else
        TI = getWindTiT(set.turb_mode, Wind.ti, collect(1:T[:nT]), SimTime)
    end

    return TI
end

