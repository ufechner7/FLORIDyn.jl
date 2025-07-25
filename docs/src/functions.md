# Exported Functions

```@meta
CurrentModule = FLORIDyn
```

# Calculating the wind directions
```@docs
getWindDirT
getWindDirT_EnKF
```

# Calculating the wind shear
```@docs
getWindShearT
```

# Calculating and applying the wind turbulence
```@docs
getWindTiT
getWindTiT_EnKF
correctTI!
```

# FLORIS
```@docs
calcCt
centerline
discretizeRotor
init_states
getVars
getPower
runFLORIS
```

# FLORIDyn
```@docs
initSimulation
findTurbineGroups
prepareSimulation
perturbationOfTheWF!
setUpTmpWFAndRun
interpolateOPs
iterateOPs!
angSOWFA2world
runFLORIDyn
```

