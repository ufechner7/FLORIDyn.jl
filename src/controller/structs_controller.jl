# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# The different controller types.

"""
    ControllerModel

An abstract type representing a controller model for wind turbines.
Subtypes of `ControllerModel` should implement specific control strategies for turbine operation.

See: [Defining the controller](@ref) for more details.
"""
abstract type ControllerModel end

"""
    Yaw_Constant <: ControllerModel

A marker struct used to represent a constant yaw control strategy.
In this mode, turbines maintain a fixed yaw angle throughout the simulation.
"""
struct Yaw_Constant <: ControllerModel end

"""
    Yaw_InterpTurbine <: ControllerModel

A marker struct used to indicate yaw control with turbine interpolation.
This mode allows for interpolated yaw angles across different turbine positions.
"""
struct Yaw_InterpTurbine <: ControllerModel end

"""
    Yaw_SOWFA <: ControllerModel

A marker struct used to represent yaw control compatible with SOWFA (Simulator fOr Wind Farm Applications).
This mode is specifically designed for integration with SOWFA simulation data.
"""
struct Yaw_SOWFA <: ControllerModel end

abstract type InductionModel end

"""
    Induction_Constant <: InductionModel

A marker struct used to represent a constant induction factor control strategy.
In this mode, all turbines maintain a fixed induction factor throughout the simulation.
The constant induction value is typically provided via the control data matrix, where
only the first element `con_induction_data[1,1]` is used for all turbines and time steps.

# Usage
This controller type is used with the `getInduction` function to provide constant
induction factors for wind farm control strategies.

# See Also
- [`getInduction`](@ref): Function for retrieving induction factors
- [`Induction_MPC`](@ref): Alternative controller for time-varying induction control
"""
struct Induction_Constant <: InductionModel end

"""
    Induction_MPC <: InductionModel

A marker struct used to represent a time-varying induction factor control strategy compatible with 
Model Predictive Control (MPC) approaches. In this mode, turbines can have different induction 
factors that vary over time according to a predefined control schedule.

The induction factors are provided via a control data matrix where:
- First column contains time values (in seconds)
- Subsequent columns contain induction factors for each turbine (dimensionless, typically 0.0-0.5)

The controller uses linear interpolation between time points with flat extrapolation for 
out-of-bounds time values.

# Usage
This controller type is used with the `getInduction` function to provide time-varying
induction factors for advanced wind farm control strategies, such as wake steering
and power optimization through axial induction control.

# Data Format
The control data matrix should have the structure:
```
[time₁  induction₁₁  induction₁₂  ...  induction₁ₙ]
[time₂  induction₂₁  induction₂₂  ...  induction₂ₙ]
[  ⋮        ⋮            ⋮        ⋱       ⋮     ]
[timeₘ  inductionₘ₁  inductionₘ₂  ...  inductionₘₙ]
```
where `m` is the number of time steps and `n` is the number of turbines.

# See Also
- [`getInduction`](@ref): Function for retrieving induction factors with interpolation
- [`Induction_Constant`](@ref): Alternative controller for constant induction control
"""
struct Induction_MPC <: InductionModel end
