```@meta
CurrentModule = FLORIDyn
```
# Simulation settings
## Introduction
The settings are defined in a `.yaml` file in the folder `data`. The function [`setup`](@ref) converts them 
into three settings structs, containing strings and numbers. The constructor [`Settings`](@ref) creates
a struct of marker types from these settings structs. 

## Abstract types
```@docs
DirModel
VelModel
```

# Settings
```@docs
setup
Settings
Settings(wind::Wind, sim::Sim)
turbineArrayProperties
getTurbineData
```