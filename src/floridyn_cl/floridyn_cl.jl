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
