# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function angSOWFA2world(deg_SOWFA)
    # Angle conversion SOWFA to world coordinates
    # deg_F = -deg_S + 270 deg
    # yaw angle defined clockwise, but for calculations counterclockwise
    deg_World = 270 - deg_SOWFA
    rad_World = deg2rad(deg_World)
    return rad_World
end

function initSimulation(wf::Union{Nothing, WindFarm}, sim::Sim)
    # Initialize the simulation or load an initialized state
    sim_init = lowercase(sim.init)

    if sim_init == "init"
        if sim.save_init_state
            # Save the initialization state to a file
            path = joinpath(sim.path_to_data, "T_init.jld2")
            @info "Saving windfield as $path ..."
            jldsave(path; wf)
        end
    elseif sim_init == "load"
        try
            data = load(joinpath(sim.path_to_data, "T_init.jld2"))
            wf = data["wf"]
        catch e
            @warn "Could not load T_init.jld2 from $(sim.path_to_data)\nWill proceed with initialized data." exception=e
        end
    end
    return wf
end

function perturbationOfTheWF!(wf, Wind)
    # perturbationOfTheWF! adds noise to the entire wind field state
    
    # Velocity
    if Wind.pertubation.vel
       wf.States_WF[:, 1] .+= Wind.pertubation.vel_sigma * randn(wf.nOP *wf.nT)
    end

    # Direction
    if Wind.pertubation.dir
       wf.States_WF[:, 2] .+= Wind.pertubation.dir_sigma * randn(wf.nOP *wf.nT)
    end

    # Turbulence Intensity
    if Wind.pertubation.ti
       wf.States_WF[:, 3] .+= Wind.pertubation.ti_sigma * randn(wf.nOP *wf.nT)
    end

    return nothing
end

function findTurbineGroups(wf, paramFLORIDyn)
    # Extract parameters from settings struct
    dw = paramFLORIDyn.deltaDW
    cw = paramFLORIDyn.deltaCW
    uw = paramFLORIDyn.deltaUW

    # Initialize outputs
    Tdep = Vector{Vector{Int64}}(undef,wf.nT)  # Equivalent of cell array in MATLAB[1][3]
    dep  = falses(wf.nT, wf.nT)
    
    R01(phi) = [cos(phi)  sin(phi); -sin(phi) cos(phi)]

    for iT in 1:wf.nT
        for iiT in 1:wf.nT
            if iiT == iT
                continue
            end

            # Closest OP from wake of iT to turbine iiT
            idx_range = (wf.StartI[iT]):(wf.StartI[iT] +wf.nOP - 1)
            distOP_iiT = sum((wf.posBase[iiT, 1:2]' .-wf.States_OP[idx_range, 1:2]).^2, dims=2)

            I = LinearIndices(distOP_iiT)
            I_op = I[argmin(distOP_iiT)]
        
            # println("I_op: ", I_op)
            I_op =wf.StartI[iT] + I_op - 1

            # Angle and relative vector
            phi = angSOWFA2world.(wf.States_WF[I_op, 2])
            r0 = wf.States_OP[I_op, 1:2] .- wf.posBase[iiT, 1:2]
            r1 = R01(phi) * r0

            # Apply dependency check
            if (-r1[1] <= uw*wf.D[iT]) && (r1[1] <= dw*wf.D[iT]) && (abs(r1[2]) <= cw*wf.D[iT])
                dep[iiT, iT] = true
            end
        end
    end

    for iT in 1:wf.nT
        # find indices where dep[iT, :] is true
        Tdep[iT] = findall(x -> x, dep[iT, :])
    end

    return Tdep
end

function interpolateOPs(wf)
    #wf.nT :: Int
    #wf.StartI :: Vector{Int}
    #wf.dep :: Vector{Vector{Int}}
    #wf.States_OP :: Matrix{Float64}
    #wf.posBase :: Matrix{Float64}
    #wf.nOP :: Int

    intOPs = Vector{Matrix{Float64}}(undef,wf.nT)  # Cell equivalent in Julia

    for iT in 1:wf.nT  # For every turbine
        intOPs[iT] = zeros(length(wf.dep[iT]), 4)

        for iiT in 1:length(wf.dep[iT])  # for every influencing turbine
            iiaT =wf.dep[iT][iiT]  # actual turbine index

            # Compute distances from OPs of turbine iiaT to current turbine
            start_idx =wf.StartI[iiaT]
            OP_positions =wf.States_OP[start_idx:(start_idx +wf.nOP - 1), 1:2]
            turb_pos =wf.posBase[iT, 1:2]

            # Euclidean distances to the turbine position
            dist = sqrt.(sum((OP_positions' .- turb_pos).^2, dims=2))
            dist = vec(dist)  # make it a flat vector

            # Indices of sorted distances
            sorted_indices = sortperm(dist)

            if sorted_indices[1] == 1
                # Closest is first OP (unlikely)
                intOPs[iT][iiT, :] = [wf.StartI[iiaT], 1.0,wf.StartI[iiaT] + 1, 0.0]
            elseif sorted_indices[1] ==wf.nOP
                # Closest is last OP (possible)
                intOPs[iT][iiT, :] = [wf.StartI[iiaT] + wf.nOP - 2, 0.0, wf.StartI[iiaT] + wf.nOP - 1, 1.0]
            else
                # Use two closest OPs for interpolation
                indOP1 =wf.StartI[iiaT] - 1 + sorted_indices[1]
                indOP2 =wf.StartI[iiaT] - 1 + sorted_indices[2]

                a =wf.States_OP[indOP1, 1:2]
                b =wf.States_OP[indOP2, 1:2]
                c =wf.posBase[iT, 1:2]

                ab = b .- a
                ac = c .- a
                d = dot(ab, ac) / dot(ab, ab)
                d = clamp(d, 0.0, 1.0)

                r1 = 1.0 - d
                r2 = d

                intOPs[iT][iiT, :] = [indOP1, r1, indOP2, r2]
            end
        end
    end

    return intOPs
end

function setUpTmpWFAndRun(set::Settings, wf, paramFLORIS, Wind)
    # Initialize outputs
    M = zeros(wf.nT, 3)
    wf.Weight = Vector{Vector{Float64}}(undef,wf.nT)
    wf.red_arr = ones(wf.nT,wf.nT)

    for iT in 1:wf.nT
        # Interpolate Wind field if needed
        iTWFState = copy(wf.States_WF[wf.StartI[iT], :])

        if hasfield(typeof(wf), :C_Vel)
            iTWFState[1] = dot(wf.C_Vel[iT, :],wf.States_WF[:, 1])
        end

        if hasfield(typeof(wf), :C_Dir)
            iTWFState[2] = dot(wf.C_Dir[iT, :],wf.States_WF[:, 2])
        end

        if isempty(wf.dep[iT])
            # Single turbine case
            T_red_arr, _, _ = runFLORIS(
                set,
                (wf.posBase[iT,:] +wf.posNac[iT,:])',
                iTWFState',
               wf.States_T[wf.StartI[iT], :]',
               wf.D[iT],
                paramFLORIS,
                Wind.shear
            )
            M[iT, :] = [T_red_arr, 0, T_red_arr *wf.States_WF[wf.StartI[iT], 1]]
           wf.red_arr[iT, iT] = T_red_arr
            continue
        end

        # Multi-turbine setup
        tmp_nT = length(wf.dep[iT]) + 1

        tmp_Tpos = repeat(wf.posBase[iT,:]' + wf.posNac[iT,:]', tmp_nT)
        tmp_WF   = repeat(iTWFState', tmp_nT)
        tmp_Tst  = repeat((wf.States_T[wf.StartI[iT], :])', tmp_nT)

        tmp_D = if wf.D[end] > 0
            vcat(wf.D[wf.dep[iT]],wf.D[iT])
        else
           wf.D
        end

        for iiT in 1:(tmp_nT - 1)
            OP1_i = Int(wf.intOPs[iT][iiT, 1])  # Index OP 1
            OP1_r = wf.intOPs[iT][iiT, 2]       # Ratio OP 1
            OP2_i = Int(wf.intOPs[iT][iiT, 3])  # Index OP 2
            OP2_r = wf.intOPs[iT][iiT, 4]       # Ratio OP 2

            OPi_l = OP1_r * wf.States_OP[OP1_i, :] + OP2_r * wf.States_OP[OP2_i, :]
            tmp_Tpos[iiT, :] = OPi_l[1:3]
            tmp_Tst[iiT, :] = OP1_r *wf.States_T[OP1_i, :] + OP2_r *wf.States_T[OP2_i, :]
            tmp_WF[iiT, :]  = OP1_r *wf.States_WF[OP1_i, :] + OP2_r *wf.States_WF[OP2_i, :]

            si = wf.StartI[wf.dep[iT][iiT]]

            if hasfield(typeof(wf), :C_Vel)
                C_weights = wf.C_Vel[iT, si:(si + wf.nOP - 1)]
                C_weights ./= sum(C_weights)
                tmp_WF[iiT, 1] = dot(C_weights, wf.States_WF[si:si + wf.nOP - 1, 1])
            end
            if hasfield(typeof(wf), :C_Dir)
                C_weights = wf.C_Dir[iT, si:(si + wf.nOP - 1)]
                C_weights ./= sum(C_weights)
                tmp_WF[iiT, 2] = dot(C_weights, wf.States_WF[si:si + wf.nOP - 1, 2])
            end

            tmp_phi = size(tmp_WF, 2) == 4 ? angSOWFA2world(tmp_WF[iiT, 4]) : angSOWFA2world(tmp_WF[iiT, 2])

            tmp_Tpos[iiT, 1] -= cos(tmp_phi) * OPi_l[4] - sin(tmp_phi) * OPi_l[5]
            tmp_Tpos[iiT, 2] -= sin(tmp_phi) * OPi_l[4] + cos(tmp_phi) * OPi_l[5]
            tmp_Tpos[iiT, 3] -= OPi_l[6]
        end

        # Run FLORIS                
        T_red_arr, T_aTI_arr, T_Ueff, T_weight = runFLORIS(set, tmp_Tpos, tmp_WF, tmp_Tst, tmp_D, paramFLORIS, Wind.shear)

        T_red = prod(T_red_arr)
        wf.red_arr[iT, vcat(wf.dep[iT], iT)] = T_red_arr
        T_addedTI = sqrt(sum(T_aTI_arr .^ 2))
        wf.Weight[iT] = T_weight

        if wf.D[end] <= 0
            dists = zeros(tmp_nT - 1)
            plot_WF = zeros(tmp_nT - 1, size(wf.States_WF, 2))
            plot_OP = zeros(tmp_nT - 1, 2)
            for iiT in 1:(tmp_nT - 1)
                OP1_i, OP1_r, OP2_i, OP2_r =wf.intOPs[iT][iiT, :]
                OPi_l = OP1_r *wf.States_OP[OP1_i, :] + OP2_r *wf.States_OP[OP2_i, :]
                plot_OP[iiT, :] = OPi_l[1:2]
                plot_WF[iiT, :] = OP1_r *wf.States_WF[OP1_i, :] + OP2_r *wf.States_WF[OP2_i, :]
                dists[iiT] = norm(OPi_l[1:2] .-wf.posBase[iT,1:2])
            end

            I = sortperm(dists)
            if length(I) == 1
                Ufree = plot_WF[I[1], 1]
                T_Ueff = T_red * Ufree
            else
                a = plot_OP[I[1], :]'
                b = plot_OP[I[2], :]'
                c =wf.posBase[iT, 1:2]'
                d = clamp((dot(b - a, c - a)) / dot(b - a, b - a), 0.0, 1.0)
                r1, r2 = 1.0 - d, d
                Ufree = r1 * plot_WF[I[1], 1] + r2 * plot_WF[I[2], 1]
                T_Ueff = T_red * Ufree
            end
        end

        M[iT, :] = [T_red, T_addedTI, T_Ueff]

        wS = sum(wf.Weight[iT])
        if wS > 0
           wf.Weight[iT] =wf.Weight[iT] ./ wS
        else
           wf.Weight[iT] .= 0.0
        end
    end

    return M, wf
end


"""
    runFLORIDyn(set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, 
                floridyn::FloriDyn, floris::Floris)

Main entry point for the FLORIDyn closed-loop simulation.

# Arguments
- `set::Settings`: Simulation settings and configuration parameters.
- `wf::WindFarm`: See: [WindFarm](@ref) simulation state, including turbine and wind farm states.
- `wind::Wind`: See: [Wind](@ref) field settings.
- `sim::Sim`: Simulation state or configuration object. See: [`Sim`](@ref)
- `con::Con`: Controller object or control parameters. See: [`Con`](@ref)
- `floridyn::FloriDyn`: Parameters specific to the FLORIDyn model. See: [`FloriDyn`](@ref)
- `floris::Floris`: Parameters specific to the FLORIS model. See: [`Floris`](@ref)

# Returns
A tuple `(wf, md, mi)` containing:
- `wf::WindFarm`: Updated simulation state with final turbine positions, wind field states, and operational point data
- `md::DataFrame`: Measurement data with columns:
  - `:Time`: Simulation time steps
  - `:ForeignReduction`: Wind speed reduction factors (%) due to wake effects from other turbines
  - `:AddedTurbulence`: Additional turbulence intensity (%) induced by upstream turbines
  - `:EffWindSpeed`: Effective wind speed (m/s) at each turbine after wake effects
  - `:FreeWindSpeed`: Free-stream wind speed (m/s) without wake interference
  - `:PowerGen`: Generated electrical power (MW) for each turbine
- `mi::Matrix`: Interaction matrix combining time data with turbine-to-turbine wake interaction coefficients 
                for each simulation step

# Description
Runs a closed-loop wind farm simulation using the FLORIDyn and FLORIS models, 
applying control strategies and updating turbine states over time.

"""
function runFLORIDyn(set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, floridyn::FloriDyn, floris::Floris)
    nT      = wf.nT
    nSim    = sim.n_sim_steps
    M       = zeros(nSim * nT, 6)
    M[:, 1] .= 1.0  # Set first column to 1
    M_int   = Vector{Matrix{Float64}}(undef, nSim)

    SimTime = sim.start_time

    for it in 1:nSim
        sim.sim_step = it

        # ========== PREDICTION ==========
        wf = iterateOPs!(set.iterate_mode, wf, sim, floris, floridyn)

        # ========== Wind Field Perturbation ==========
        perturbationOfTheWF!(wf, wind)

        # ========== Get FLORIS reductions ==========
        wf.dep = findTurbineGroups(wf, floridyn)
        wf.intOPs = interpolateOPs(wf)
        a, b = setUpTmpWFAndRun(set, wf, floris, wind)
        tmpM, wf = a, b
        M[(it-1)*nT+1 : it*nT, 2:4] .= tmpM
        M[(it-1)*nT+1 : it*nT, 1]   .= SimTime
        wf.States_T[wf.StartI, 3] = tmpM[:, 2]
        M_int[it] = wf.red_arr

        # ========== wind field corrections ==========
        wf, wind = correctVel(set.cor_vel_mode, set, wf, wind, SimTime, floris, tmpM)
        correctDir!(set.cor_dir_mode, set, wf, wind, SimTime)
        wf = correctTi(set.cor_turb_mode, set, wf, wind, SimTime)

        # Save free wind speed as measurement
        M[(it-1)*nT+1 : it*nT, 5] = wf.States_WF[wf.StartI, 1]

        # ========== Get Control settings ==========
        wf.States_T[wf.StartI, 2] = (
            wf.States_WF[wf.StartI, 2] .-
                getYaw(set.control_mode, con.yaw_data, collect(1:nT), SimTime)'
        )

        # ========== Calculate Power ==========
        P = getPower(wf, tmpM, floris, con)
        M[(it-1)*nT+1:it*nT, 6] = P

        SimTime += sim.time_step
    end
    # Convert `M` to DataFrame and scale measurements
    md = DataFrame(
        (M * diagm([1; 100; 100; 1; 1; 1e-6])),
        [:Time, :ForeignReduction, :AddedTurbulence, :EffWindSpeed, :FreeWindSpeed, :PowerGen]
    )
    mi = hcat(md.Time, hcat(M_int...)')
    return wf, md, mi
end
