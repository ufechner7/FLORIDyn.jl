# Exported Types

```@meta
CurrentModule = FLORIDyn
```

## Markers for defining the wind direction
An instance of these structs needs to be passed to the functions that calculate the wind direction.

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

## Markers for defining the wind shear
An instance of these structs needs to be passed to the functions that calculate the wind shear.
```@docs
Shear_Interpolation
Shear_PowerLaw
```

## Markers for defining the wind turbulence
An instance of these structs needs to be passed to the functions that calculate the wind turbulence.
```@docs
 TI_Constant
 TI_EnKF_InterpTurbine
 TI_Interpolation
 TI_InterpTurbine
```