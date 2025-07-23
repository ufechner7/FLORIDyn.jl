# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function correctTi(::TI_None, set, wf, Wind, SimTime)
    # correctTi updates the turbulent intensity (TI) value inT.States_WF
    # at the rowT.StartI and column 3, using data from getDataTI.

    # Get new turbulent intensity value
    TI = getDataTI(set, Wind, wf, SimTime)

    # Update the TI value in the States_WF matrix
    try
       wf.States_WF[wf.StartI, 3] .= TI'
    catch e
        error("Error updatingT.States_WF: $(e.msg)")
    end

    return wf
end

function getDataTI(set, Wind, wf, SimTime)
    # GETDATATI retrieves the data for the ambient turbulence intensity

    if Wind.input_ti == "EnKF_InterpTurbine"
        TI = getWindTiT_EnKF(set.turb_mode, Wind.ti, collect(1:wf.nT), SimTime)
    elseif Wind.input_ti == "EnKF_ZOH"
        TI =wf.States_WF[wf.StartI, 3]
    elseif Wind.input_ti == "EnKF_RW"
        TI =wf.States_WF[wf.StartI, 3]
        # 'phi' is not defined in your snippet. Replace with wf.nT if that's correct.
        TI = TI .+ (randn(1, wf.nT) * Wind.ti.CholSig)'
    elseif Wind.input_ti == "CLC_weighted_ZOH"
        TI = wf.C_TI *wf.States_WF[:, 3]
    else
        TI = getWindTiT(set.turb_mode, Wind.ti, collect(1:wf.nT), SimTime)
    end

    return TI
end

