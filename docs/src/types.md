# Exported Types

```@meta
CurrentModule = FLORIDyn
```

## Marker types for defining the wind direction
An instance of these structs needs to be passed to the functions that calculate the wind direction.
Example:
```julia
dir_mode = Direction_constant()
```
```@docs
Direction_Constant
Direction_Constant_wErrCov
```