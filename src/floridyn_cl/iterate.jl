# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function iterateOPs!(::IterateOPs_basic, T, Sim, paramFLORIS, paramFLORIDyn)
    # Save turbine OPs
    tmpOPStates = copy(T[:States_OP][T[:StartI], :])
    tmpTStates  = copy(T[:States_T][T[:StartI], :])
    tmpWFSTates = copy(T[:States_WF][T[:StartI], :])

    # Shift states
    # Downwind step
    step_dw = Sim.time_step .*T.States_WF[:, 1] .* Sim.dyn.advection
   T.States_OP[:, 4] .+= step_dw

    # Crosswind step
    deflection = centerline(T[:States_OP],T.States_T,T.States_WF, paramFLORIS,T.D[1])
    step_cw = deflection .-T.States_OP[:, 5:6]
   T.States_OP[:, 5:6] .= deflection

    # World coordinate system adjustment
    phiW = angSOWFA2world.(T[:States_WF][:, 2])
   T.States_OP[:, 1] .+= cos.(phiW) .* step_dw .- sin.(phiW) .* step_cw[:, 1]
   T.States_OP[:, 2] .+= sin.(phiW) .* step_dw .+ cos.(phiW) .* step_cw[:, 1]
   T.States_OP[:, 3] .+= step_cw[:, 2]

    # Circshift & init first OPs
    # OPs
   T.States_OP = circshift(T[:States_OP], (1, 0))
    #T.States_OP = circshift(T[:States_OP], (1, 0))
   T.States_OP[T[:StartI], :] = tmpOPStates

    # Turbines
   T.States_T = circshift(T[:States_T], (1, 0))
   T.States_T[T[:StartI], :] = tmpTStates

    # Wind Farm
   T.States_WF = circshift(T[:States_WF], (1, 0))
   T.States_WF[T[:StartI], :] = tmpWFSTates

    # Check if OPs are in order
    for iT in 1:T[:nT]
        inds =T.StartI[iT]:(T[:StartI][iT] +T.nOP - 1)

        indOP = sortperm(T[:States_OP][inds, 4])
        if indOP != sort(indOP)  # check if already sorted
           T.States_OP[inds, :] =T.States_OP[inds[indOP], :]
           T.States_T[inds, :]  =T.States_T[inds[indOP], :]
           T.States_WF[inds, :] =T.States_WF[inds[indOP], :]
        end

    end
    return T
end
