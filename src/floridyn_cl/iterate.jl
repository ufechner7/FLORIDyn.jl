# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function iterateOPs!(::IterateOPs_basic, wf, Sim, paramFLORIS, paramFLORIDyn)
    # Save turbine OPs
    tmpOPStates = copy(wf.States_OP[wf.StartI, :])
    tmpTStates  = copy(wf.States_T[wf.StartI, :])
    tmpWFSTates = copy(wf.States_WF[wf.StartI, :])

    # Shift states
    # Downwind step
    step_dw = Sim.time_step .* wf.States_WF[:, 1] .* Sim.dyn.advection
    wf.States_OP[:, 4] .+= step_dw

    # Crosswind step
    deflection = centerline(wf.States_OP, wf.States_T, wf.States_WF, paramFLORIS, wf.D[1])
    step_cw = deflection .- wf.States_OP[:, 5:6]
    wf.States_OP[:, 5:6] .= deflection

    # World coordinate system adjustment
    phiW = angSOWFA2world.(wf.States_WF[:, 2])
    wf.States_OP[:, 1] .+= cos.(phiW) .* step_dw .- sin.(phiW) .* step_cw[:, 1]
    wf.States_OP[:, 2] .+= sin.(phiW) .* step_dw .+ cos.(phiW) .* step_cw[:, 1]
    wf.States_OP[:, 3] .+= step_cw[:, 2]

    # Circshift & init first OPs
    # OPs
    wf.States_OP = circshift(wf.States_OP, (1, 0))
    wf.States_OP[wf.StartI, :] = tmpOPStates

    # Turbines
    wf.States_T = circshift(wf.States_T, (1, 0))
    wf.States_T[wf.StartI, :] = tmpTStates

    # Wind Farm
    wf.States_WF = circshift(wf.States_WF, (1, 0))
    wf.States_WF[wf.StartI, :] = tmpWFSTates

    # Check if OPs are in order
    for iT in 1:wf.nT
        inds =wf.StartI[iT]:(wf.StartI[iT] +wf.nOP - 1)

        indOP = sortperm(wf.States_OP[inds, 4])
        if indOP != sort(indOP)  # check if already sorted
           wf.States_OP[inds, :] =wf.States_OP[inds[indOP], :]
           wf.States_T[inds, :]  =wf.States_T[inds[indOP], :]
           wf.States_WF[inds, :] =wf.States_WF[inds[indOP], :]
        end

    end
    return wf
end
