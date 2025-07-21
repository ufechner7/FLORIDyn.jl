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
    #     T.intOPs = interpolateOPs(T)
    #     tmpM, T = setUpTmpWFAndRun(T, paramFLORIS, Wind)
    #     M[(it-1)*nT+1 : it*nT, 2:4] = tmpM
    #     M[(it-1)*nT+1 : it*nT, 1] = SimTime
    #     T.States_T[T[:StartI], 3] = tmpM[:, 2]
    #     M_int[it] = T.red_arr

    #     # ========== Wind field corrections ==========
    #     T, Wind = correctVel(T, Wind, SimTime, paramFLORIS, tmpM)
    #     T = correctDir(T, Wind, SimTime)
    #     T = correctTi(T, Wind, SimTime)

    #     # Save free wind speed as measurement
    #     M[(it-1)*nT+1 : it*nT, 5] = T[:States_WF][T[:StartI], 1]

    #     # ========== Get Control settings ==========
    #     T.States_T[T[:StartI], 2] = (
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
