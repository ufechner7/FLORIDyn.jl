# High-Level Functions

```@meta
CurrentModule = FLORIDyn
```

# FLORIS
```@docs
calcCt
centerline!
discretizeRotor
init_states
getVars!
getPower
runFLORIS!
```

Note: `runFLORIS!` is allocation-free and returns nothing. Read results from the
`FLORISBuffers` you passed in (fields `T_red_arr`, `T_aTI_arr`, `T_Ueff`, `T_weight`).

# FLORIDyn
```@docs
initSimulation
findTurbineGroups
prepareSimulation
perturbationOfTheWF!
setUpTmpWFAndRun!
interpolateOPs!
iterateOPs!
angSOWFA2world
runFLORIDyn
run_floridyn
```

# Visualization
```@docs
getMeasurements
create_thread_buffers
getMeasurementsP
calcFlowField
plotFlowField
plot_flow_field
plotMeasurements
plot_measurements
close_all
```

# Video Creation
```@docs
cleanup_video_folder
createVideo
createAllVideos
natural_sort_key
```

