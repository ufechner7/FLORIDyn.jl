# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# the different wind turbulence types (turb_mode)
"""
    TI_Constant

A marker struct representing a constant turbulence intensity. 
"""
struct TI_Constant end

"""
    TI_EnKF_InterpTurbine

A marker struct representing the Turbulence Intensity (TI) Ensemble Kalman Filter (EnKF) interpolation model.
"""
struct TI_EnKF_InterpTurbine end

"""
    TI_Interpolation

A marker struct representing the interpolation method for modeling the turbulence.
"""
struct TI_Interpolation end


"""
    TI_InterpTurbine

A marker struct representing an interpolated turbine model for turbulence intensity calculations.
"""
struct TI_InterpTurbine end
