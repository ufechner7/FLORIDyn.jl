# Exported Types

```@meta
CurrentModule = FLORIDyn
```
## The Wind Farm Simulation struct
```@docs
WindFarm
```
## Abstract types
```@docs
VelModel
DirModel
ShearModel
TurbulenceModel
VelCorrection
DirCorrection
TurbulenceCorrection
IterateOPs_model
ControllerModel
```

## Defining the wind velocity model
An instance of these structs needs to be passed to the functions that calculate the wind velocity. They are
all subtypes of [`VelModel`](@ref)
```@docs
Velocity_Constant
Velocity_Constant_wErrorCov
Velocity_EnKF_InterpTurbine
Velocity_I_and_I
Velocity_Interpolation
Velocity_Interpolation_wErrorCov
Velocity_InterpTurbine
Velocity_InterpTurbine_wErrorCov
Velocity_RW_with_Mean
Velocity_ZOH_wErrorCov
```

## Defining the wind direction model
An instance of these structs needs to be passed to the functions that calculate the wind direction. They are
all subtypes of [`DirModel`](@ref)

```@docs
Direction_Constant
Direction_Constant_wErrorCov
Direction_EnKF_InterpTurbine
Direction_Interpolation
Direction_Interpolation_wErrorCov
Direction_InterpTurbine
Direction_InterpTurbine_wErrorCov
Direction_RW_with_Mean
```

## Defining the wind shear model
An instance of these structs needs to be passed to the functions that calculate the wind shear.
```@docs
Shear_Interpolation
Shear_PowerLaw
```

## Defining the wind turbulence model
An instance of these structs needs to be passed to the functions that calculate the wind turbulence.
```@docs
 TI_Constant
 TI_EnKF_InterpTurbine
 TI_Interpolation
 TI_InterpTurbine
```

## Defining the velocity correction
An instance of these structs needs to be passed to the functions that correct wind velocity. They are
all subtypes of [`VelCorrection`](@ref)

```@docs
Velocity_Influence
Velocity_None
Velocity_wGaspariAndCohn
```

## Defining the direction correction
An instance of these structs needs to be passed to the functions that correct wind direction. They are
all subtypes of [`DirCorrection`](@ref)

```@docs
Direction_All
Direction_Influence
Direction_None
```

## Defining the turbulence correction
An instance of these structs needs to be passed to the functions that correct turbulence intensity. They are
all subtypes of [`TurbulenceCorrection`](@ref)

```@docs
TI_Influence
TI_None
TI_wGaspariAndCohn
```

## Defining the controller
An instance of these structs needs to be passed to the functions that control turbine behavior. They are all subtypes of [`ControllerModel`](@ref)

```@docs
Yaw_Constant
Yaw_InterpTurbine
Yaw_SOWFA
```

## Defining the OP iteration model
An instance of these structs needs to be passed to the functions that iterate operational points (OPs) through the wind field. They are all subtypes of [`IterateOPs_model`](@ref)

```@docs
IterateOPs_average
IterateOPs_basic
IterateOPs_buffer
IterateOPs_maximum
IterateOPs_weighted
```

## Types for storing wind field params
```@docs
WindDirType
WindDirMatrix
WindDirTriple
WindVelType
WindVelMatrix
WindShear
WindPerturbation
```

## Other types
```@docs
States
IterateOPsBuffers
```
