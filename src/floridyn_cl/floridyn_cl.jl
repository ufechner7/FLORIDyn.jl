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

function initSimulation(T, Wind, Sim, Con, paramFLORIDyn, paramFLORIS)
    # Initialize the simulation or load an initialized state
    sim_init = lowercase(Sim.init)

    if sim_init == "init"
        if getfield(Sim, :SaveInitState) # or Sim.SaveInitState if Sim is mutable struct
            # Save the initialization state to a file
            jldsave(joinpath(Sim.path_to_data, "T_init.jld2"); T)
        end
    elseif sim_init == "load"
        try
            data = load(joinpath(Sim.path_to_data, "T_init.jld2"))
            T = data["T"]
        catch e
            @warn "Could not load T_init.jld2 from $(Sim.PathToSim)\nWill proceed with initialized data." exception=e
        end
    end
    return T
end

function pertubationOfTheWF!(T, Wind)
    # pertubationOfTheWF! adds noise to the entire wind field state
    
    # Velocity
    if Wind.pertubation.vel
        T[:States_WF][:, 1] .+= Wind.pertubation.vel_sigma * randn(T[:nOP] * T[:nT])
    end

    # Direction
    if Wind.pertubation.dir
        T[:States_WF][:, 2] .+= Wind.pertubation.dir_sigma * randn(T[:nOP] * T[:nT])
    end

    # Turbulence Intensity
    if Wind.pertubation.ti
        T[:States_WF][:, 3] .+= Wind.pertubation.ti_sigma * randn(T[:nOP] * T[:nT])
    end

    return nothing
end

function findTurbineGroups(T, paramFLORIDyn)
    # Extract parameters from settings struct
    dw = paramFLORIDyn.deltaDW
    cw = paramFLORIDyn.deltaCW
    uw = paramFLORIDyn.deltaUW

    # Initialize outputs
    Tdep = Vector{Any}(undef, T[:nT])  # Equivalent of cell array in MATLAB[1][3]
    dep = falses(T[:nT], T[:nT])
    
    R01(phi) = [cos(phi)  sin(phi); -sin(phi) cos(phi)]

    for iT in 1:T[:nT]
        for iiT in 1:T[:nT]
            if iiT == iT
                continue
            end

            # Closest OP from wake of iT to turbine iiT
            idx_range = (T[:StartI][iT]):(T[:StartI][iT] + T[:nOP] - 1)
            distOP_iiT = sum((T[:posBase][iiT, 1:2]' .- T[:States_OP][idx_range, 1:2]).^2, dims=2)

            I = LinearIndices(distOP_iiT)
            I_op = I[argmin(distOP_iiT)]
        
            # println("I_op: ", I_op)
            I_op = T[:StartI][iT] + I_op - 1

            # Angle and relative vector
            phi = angSOWFA2world.(T[:States_WF][I_op, 2])
            r0 = T[:States_OP][I_op, 1:2] .- T[:posBase][iiT, 1:2]
            r1 = R01(phi) * r0

            # Apply dependency check
            if (-r1[1] <= uw*T[:D][iT]) && (r1[1] <= dw*T[:D][iT]) && (abs(r1[2]) <= cw*T[:D][iT])
                dep[iiT, iT] = true
            end
        end
    end

    for iT in 1:T[:nT]
        # find indices where dep[iT, :] is true
        Tdep[iT] = findall(x -> x, dep[iT, :])
    end

    return Tdep
end

function interpolateOPs(T)
    # T[:nT] :: Int
    # T[:StartI] :: Vector{Int}
    # T[:dep] :: Vector{Vector{Int}}
    # T[:States_OP] :: Matrix{Float64}
    # T[:posBase] :: Matrix{Float64}
    # T[:nOP] :: Int

    intOPs = Vector{Matrix{Float64}}(undef, T[:nT])  # Cell equivalent in Julia

    for iT in 1:T[:nT]  # For every turbine
        intOPs[iT] = zeros(length(T[:dep][iT]), 4)

        for iiT in 1:length(T[:dep][iT])  # for every influencing turbine
            iiaT = T[:dep][iT][iiT]  # actual turbine index

            # Compute distances from OPs of turbine iiaT to current turbine
            start_idx = T[:StartI][iiaT]
            OP_positions = T[:States_OP][start_idx:start_idx + T[:nOP] - 1, 1:2]
            turb_pos = T[:posBase][iT, 1:2]

            # Euclidean distances to the turbine position
            dist = sqrt.(sum((OP_positions' .- turb_pos).^2, dims=2))
            dist = vec(dist)  # make it a flat vector

            # Indices of sorted distances
            sorted_indices = sortperm(dist)

            if sorted_indices[1] == 1
                # Closest is first OP (unlikely)
                intOPs[iT][iiT, :] = [T[:StartI][iiaT], 1.0, T[:StartI][iiaT] + 1, 0.0]
            elseif sorted_indices[1] == T[:nOP]
                # Closest is last OP (possible)
                intOPs[iT][iiT, :] = [T[:StartI][iiaT] + T[:nOP] - 2, 0.0, T[:StartI][iiaT] + T[:nOP] - 1, 1.0]
            else
                # Use two closest OPs for interpolation
                indOP1 = T[:StartI][iiaT] - 1 + sorted_indices[1]
                indOP2 = T[:StartI][iiaT] - 1 + sorted_indices[2]

                a = T[:States_OP][indOP1, 1:2]
                b = T[:States_OP][indOP2, 1:2]
                c = T[:posBase][iT, 1:2]

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

function setUpTmpWFAndRun(set::Settings, T, paramFLORIS, Wind)
    # Initialize outputs
    M = zeros(T[:nT], 3)
    T[:Weight] = Vector{Any}(undef, T[:nT])
    T[:red_arr] = ones(T[:nT], T[:nT])

    for iT in 1:T[:nT]
        # Interpolate Wind field if needed
        iTWFState = copy(T[:States_WF][T[:StartI][iT], :])

        if hasfield(typeof(T), :C_Vel)
            iTWFState[1] = dot(T[:C_Vel][iT, :], T[:States_WF][:, 1])
        end

        if hasfield(typeof(T), :C_Dir)
            iTWFState[2] = dot(T.C_Dir[iT, :], T[:States_WF][:, 2])
        end

        if isempty(T[:dep][iT])
            # Single turbine case
            T_red_arr, _, _ = runFLORIS(
                set,
                T[:posBase][iT,:] + T[:posNac][iT,:],
                iTWFState,
                T[:States_T][T[:StartI][iT], :],
                T[:D][iT],
                paramFLORIS,
                Wind.shear
            )
            M[iT, :] = [T_red_arr, 0, T_red_arr * T[:States_WF][T[:StartI][iT], 1]]
            T[:red_arr][iT, iT] = T_red_arr
            continue
        end

        # Multi-turbine setup
        tmp_nT = length(T[:dep][iT]) + 1

        tmp_Tpos = repeat(T[:posBase][iT,:]' + T[:posNac][iT,:]', tmp_nT)
        tmp_WF   = repeat(iTWFState', tmp_nT)
        tmp_Tst  = repeat((T[:States_T][T[:StartI][iT], :])', tmp_nT)

        tmp_D = if T[:D][end] > 0
            vcat(T[:D][T[:dep][iT]], T[:D][iT])
        else
            T[:D]
        end

        for iiT in 1:(tmp_nT - 1)
            OP1_i = Int(T[:intOPs][iT][iiT, 1])  # Index OP 1
            OP1_r = T[:intOPs][iT][iiT, 2]       # Ratio OP 1
            OP2_i = Int(T[:intOPs][iT][iiT, 3])  # Index OP 2
            OP2_r = T[:intOPs][iT][iiT, 4]       # Ratio OP 2
            # println("OP1_i: ", OP1_i, " OP2_i: ", OP2_i, " iiT: ", iiT)
            # println("OP1_r: ", OP1_r, " OP2_r: ", OP2_r)
            # OP1_i, OP1_r, OP2_i, OP2_r = T[:intOPs][iT][iiT, :]  # Assumes row-major

            OPi_l = OP1_r * T[:States_OP][OP1_i, :] + OP2_r * T[:States_OP][OP2_i, :]
            println("OPi_l: ", OPi_l)
            println("tmp_Tpos[iiT, :]: ", tmp_Tpos[iiT, :])
            tmp_Tpos[iiT, :] = OPi_l[1:3]
            tmp_Tst[iiT, :] = OP1_r * T[:States_T][OP1_i, :] + OP2_r * T[:States_T][OP2_i, :]
            tmp_WF[iiT, :]  = OP1_r * T[:States_WF][OP1_i, :] + OP2_r * T[:States_WF][OP2_i, :]

            si = T[:StartI][T[:dep][iT][iiT]]

            if hasfield(typeof(T), :C_Vel)
                C_weights = T[:C_Vel][iT, si:(si + T[:nOP] - 1)]
                C_weights ./= sum(C_weights)
                tmp_WF[iiT, 1] = dot(C_weights, T[:States_WF][si:si + T[:nOP] - 1, 1])
            end
            if hasfield(typeof(T), :C_Dir)
                C_weights = T.C_Dir[iT, si:(si + T[:nOP] - 1)]
                C_weights ./= sum(C_weights)
                tmp_WF[iiT, 2] = dot(C_weights, T[:States_WF][si:si + T[:nOP] - 1, 2])
            end

            tmp_phi = size(tmp_WF, 2) == 4 ? angSOWFA2world(tmp_WF[iiT, 4]) : angSOWFA2world(tmp_WF[iiT, 2])

            tmp_Tpos[iiT, 1] -= cos(tmp_phi) * OPi_l[4] - sin(tmp_phi) * OPi_l[5]
            tmp_Tpos[iiT, 2] -= sin(tmp_phi) * OPi_l[4] + cos(tmp_phi) * OPi_l[5]
            tmp_Tpos[iiT, 3] -= OPi_l[6]
        end

        # Run FLORIS
        T_red_arr, T_aTI_arr, T_Ueff, T_weight = runFLORIS(set, tmp_Tpos, tmp_WF, tmp_Tst, tmp_D, paramFLORIS, Wind.shear)

        T_red = prod(T_red_arr)
        T[:red_arr][iT, vcat(T[:dep][iT], iT)] = T_red_arr
        T_addedTI = sqrt(sum(T_aTI_arr .^ 2))
        T[:Weight][iT] = T_weight

        if T[:D][end] <= 0
            dists = zeros(tmp_nT - 1)
            plot_WF = zeros(tmp_nT - 1, size(T[:States_WF], 2))
            plot_OP = zeros(tmp_nT - 1, 2)
            for iiT in 1:(tmp_nT - 1)
                OP1_i, OP1_r, OP2_i, OP2_r = T[:intOPs][iT][iiT, :]
                OPi_l = OP1_r * T[:States_OP][OP1_i, :] + OP2_r * T[:States_OP][OP2_i, :]
                plot_OP[iiT, :] = OPi_l[1:2]
                plot_WF[iiT, :] = OP1_r * T[:States_WF][OP1_i, :] + OP2_r * T[:States_WF][OP2_i, :]
                dists[iiT] = norm(OPi_l[1:2] .- T[:posBase][iT,1:2])
            end

            I = sortperm(dists)
            if length(I) == 1
                Ufree = plot_WF[I[1], 1]
                T_Ueff = T_red * Ufree
            else
                a = plot_OP[I[1], :]'
                b = plot_OP[I[2], :]'
                c = T[:posBase][iT, 1:2]'
                d = clamp((dot(b - a, c - a)) / dot(b - a, b - a), 0.0, 1.0)
                r1, r2 = 1.0 - d, d
                Ufree = r1 * plot_WF[I[1], 1] + r2 * plot_WF[I[2], 1]
                T_Ueff = T_red * Ufree
            end
        end

        M[iT, :] = [T_red, T_addedTI, T_Ueff]

        wS = sum(T[:Weight][iT])
        if wS > 0
            T[:Weight][iT] = T[:Weight][iT] ./ wS
        else
            T[:Weight][iT] .= 0.0
        end
    end

    return M, T
end


function FLORIDynCL(set::Settings, T, Wind, Sim, Con, paramFLORIDyn, paramFLORIS)
    # OUTPUTS:
    # T := Simulation state (OP states, Turbine states, wind field states(OPs))
    # Mt := Measurements from the simulation (Power, tbd)
    # Mint := Interaction matrices for the turbines

    nT      = T[:nT]
    nSim    = Sim.n_sim_steps
    M       = zeros(nSim * nT, 6)
    M[:, 1] .= 1.0  # Set first column to 1
    M_int   = Vector{Any}(undef, nSim)

    SimTime = Sim.start_time

    for it in 1:nSim
        Sim.sim_step = it

        # ========== PREDICTION ==========
        T = iterateOPs!(set.iterate_mode, T, Sim, paramFLORIS, paramFLORIDyn)

        # ========== Wind Field Perturbation ==========
        pertubationOfTheWF!(T, Wind)

        # ========== Get FLORIS reductions ==========
        T[:dep] = findTurbineGroups(T, paramFLORIDyn)
        T[:intOPs] = interpolateOPs(T)
        tmpM, T = setUpTmpWFAndRun(set, T, paramFLORIS, Wind)
    #     M[(it-1)*nT+1 : it*nT, 2:4] = tmpM
    #     M[(it-1)*nT+1 : it*nT, 1] = SimTime
    #     T[:States_T][T[:StartI], 3] = tmpM[:, 2]
    #     M_int[it] = T[:red_arr]

    #     # ========== Wind field corrections ==========
    #     T, Wind = correctVel(T, Wind, SimTime, paramFLORIS, tmpM)
    #     T = correctDir(T, Wind, SimTime)
    #     T = correctTi(T, Wind, SimTime)

    #     # Save free wind speed as measurement
    #     M[(it-1)*nT+1 : it*nT, 5] = T[:States_WF][T[:StartI], 1]

    #     # ========== Get Control settings ==========
    #     T[:States_T][T[:StartI], 2] = (
    #         T[:States_WF][T[:StartI], 2] .-
    #         getYaw(Con.YawData, collect(1:nT), SimTime)'
    #     )

    #     # ========== Calculate Power ==========
    #     P = getPower(T, tmpM, paramFLORIS, Con)
    #     M[(it-1)*nT+1:it*nT, 6] = P

        SimTime += Sim.time_step
    end
    # # Convert `M` to DataFrame and scale measurements
    # Mt = DataFrame(
    #     (M .* [1; 100; 100; 1; 1; 1e-6])',
    #     [:Time, :ForeignReduction, :AddedTurbulence, :EffWindSpeed, :FreeWindSpeed, :PowerGen]
    # )
    # Mint = hcat(Mt.Time, hcat(M_int...)')
    Mt = nothing
    Mint = nothing
    return T, Mt, Mint
end
