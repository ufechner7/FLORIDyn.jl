# Exported Functions

```@meta
CurrentModule = FLORIDyn
```

# Calculating the wind directions
```@docs
getWindDirT
getWindDirT_EnKF
correctDir!
```

# Calculating the wind shear
```@docs
getWindShearT
```

# Calculating the wind turbulence
```@docs
getDataTI
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

