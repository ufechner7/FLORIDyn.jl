# Exported Types

```@meta
CurrentModule = FLORIDyn
```
## The Wind Farm Simulation struct
```@docs
WindFarm
```

## Defining the wind velocity model
An instance of these structs needs to be passed to the functions that calculate the wind velocity. They are
all subtypes of [`VelModel`](@ref)
```@docs
Velocity_Constant
Velocity_Constant_wErrorCov
Velocity_EnKF_InterpTurbine
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
Direction_wGaspariAndCohn
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
```
