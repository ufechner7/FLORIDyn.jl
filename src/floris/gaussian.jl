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
