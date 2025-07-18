# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    CalcCt(a, _)

Calculate the thrust coefficient (Ct) for a wind turbine based on the axial induction factor `a`.

# Arguments
- `a::Number`: Axial induction factor, typically between 0 and 0.5.
- _: unused parameter

# Returns
- `Ct::Number`: The calculated thrust coefficient.
"""
function CalcCt(a, _)
    Ct = 4 .* a .* (1 .- a)
    return Ct
end

mutable struct States
    T_names::Vector{String}
    Turbine::Int
    OP_names::Vector{String}
    OP::Int
    WF_names::Vector{String}
    WF::Int
end

function States()
    # Turbine states
    T_names = ["a", "yaw", "TI"]
    Turbine = length(T_names)

    # Observation point states
    OP_names = ["x0", "y0", "z0", "x1", "y1", "z1"]
    OP = length(OP_names)

    # Wind field states
    WF_names = ["wind_vel", "wind_dir", "TI0"]
    WF = length(WF_names)

    return States(T_names, Turbine, OP_names, OP, WF_names, WF)
end

