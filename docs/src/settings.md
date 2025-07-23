```@meta
CurrentModule = FLORIDyn
```
# Simulation settings
## Introduction
The settings are defined in a `.yaml` file in the folder `data`. The function [`setup`](@ref) converts them into three settings structs, containing strings and numbers. The constructor [`Settings`](@ref) creates a struct of marker types from these settings structs. 

For each of these yaml files there must be a folder with the same name (without the `.yaml` extension) in the data directory. In this folder the required `.csv` files must be stored.

## Abstract types
```@docs
VelModel
DirModel
ShearModel
TurbulenceModel
VelCorrection
DirCorrection
ControllerModel
```

# Types created from the yaml file
```@docs
Sim
Wind
Con
Floris
FloriDyn
```

# Settings
```@docs
setup
Settings
Settings(wind::Wind, sim::Sim, con::Con)
turbineArrayProperties
getTurbineData
```