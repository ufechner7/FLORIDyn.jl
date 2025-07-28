# Exported Functions

```@meta
CurrentModule = FLORIDyn
```

# Calculating the wind directions
```@docs
getWindDirT
getWindDirT_EnKF
getDataDir
correctDir!
```

# Calculating the wind velocity
```@docs
getDataVel
```

# Calculating the wind shear
```@docs
getWindShearT
```

# Calculating the wind turbulence
```@docs
getDataTI
getWindTiT
correctTI!
```

# Controller functions
```@docs
getYaw
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

# Visualization
```@docs
calcFlowField
plotFlowField
```

