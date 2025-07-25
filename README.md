# FLORIDyn

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ufechner7.github.io/FLORIDyn.jl/dev)
[![Build Status](https://github.com/ufechner7/FLORIDyn.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ufechner7/FLORIDyn.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/ufechner7/FLORIDyn.jl/graph/badge.svg?token=O7wXT62VSR)](https://codecov.io/gh/ufechner7/FLORIDyn.jl)

## Introduction
A dynamic wind farm simulation software, translated from https://github.com/TUDelft-DataDrivenControl/FLORIDyn_Matlab, which was written by Marcus Becker.

His code uses one engineering model from the quasi static wind farm simulation software [FLORIS](https://github.com/NREL/floris), developed by NREL.

## Model features
- Simulate wind farms dynamically at a low computational cost
- Estimate the power generated, added turbulence, and wake-induced losses.
- Apply heterogeneous and time-varying wind speeds and directions
- Test different modeling approaches

## Status:
The basic example works.
```julia
include("examples/main.jl")
```
### TODO
- implement the visualization
- add more unit tests
- implement the missing functions

A Python version of FLORIDyn is available at https://github.com/TUDelft-DataDrivenControl/OFF .

## Installation
Install [Julia 1.10](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html) or later, if you haven't already.

Before installing this software it is suggested to create a new project, for example like this:
```bash
mkdir test
cd test
julia --project="."
```
Then add FLORIDyn from  Julia's package manager, by typing:
```julia
using Pkg
pkg"add https://github.com/ufechner7/FLORIDyn.jl"
``` 
at the Julia prompt. You can run the unit tests with the command:
```julia
pkg"test FLORIDyn"
```

## Documentation
The documentation is available [here](https://ufechner7.github.io/FLORIDyn.jl/dev/).

## References
Citation of the FLORIDyn model:
FLORIDyn - A dynamic and flexible framework for real-time wind farm control, M. Becker, D. Allaerts, J.W. van Wingerden, 2022, http://doi.org/10.1088/1742-6596/2265/3/032103

Used FLORIS model:
Experimental and theoretical study of wind turbine wakes in yawed conditions, M. Bastankhah, F. Porté-Agel, 2020, http://doi.org/10.1017/jfm.2016.595

Additional references for smaller subcomponents can be found in the code or in the related publications.

## Publications that use FLORIDyn
- Ensemble-Based Flow Field Estimation Using the Dynamic Wind Farm Model FLORIDyn, M. Becker, D. Allaerts, J.W. van Wingerden, 2022, http://doi.org/10.3390/en15228589
- Wind pattern clustering of high frequent field measurements for dynamic wind farm flow control, M. Becker, D. Allaerts, J.W. van Wingerden, 2024, http://doi.org/10.1088/1742-6596/2767/3/032028 
- Sensitivity analysis and Bayesian calibration of a dynamic wind farm control model: FLORIDyn, V.V. Dighe, M. Becker, wf. Göçmen, B. Sanderse, J.W. van Wingerden, 2022, http://doi.org/10.1088/1742-6596/2265/2/022062
- Time-shifted cost function design for more efficient dynamic wind farm flow control, M. Becker, D. Allaerts, J.W. van Wingerden, 2024, http://doi.org/10.1109/CCTA60707.2024.10666535
- Suitability of Dynamic Wake Models for AEP Estimation: A Wind Farm-Scale Validation Study, M. Van der Straeten, http://resolver.tudelft.nl/uuid:f35617a2-2409-439b-8bc2-6334b807ce1f 
- Scaling DMD modes for modeling Dynamic Induction Control wakes in various wind speeds, J. Gutknecht, M. Becker, C. Muscari, wf. Lutz, J.W. van Wingerden, 2023, http://doi.org/10.1109/CCTA54093.2023.10252400
- Model predictive control of wakes for wind farm power tracking, A. Sterle, C.A. Hans, J. Raisch, 2024, http://doi.org/10.1088/1742-6596/2767/3/032005

