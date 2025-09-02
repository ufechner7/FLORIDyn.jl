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
prepare_rotor_points!
handle_single_turbine!
setup_computation_buffers!
compute_wake_effects!
compute_final_wind_shear!
```

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
create_thread_buffers
getMeasurements
calcFlowField
calc_rel_power
plotFlowField
plot_flow_field
plotMeasurements
plot_measurements
plot_x
prepare_large_plot_inputs
close_all
```

# Video Creation
```@docs
cleanup_video_folder
createVideo
createAllVideos
natural_sort_key
```

