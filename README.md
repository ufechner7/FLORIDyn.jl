# FLORIDyn

[![Build Status](https://github.com/ufechner7/FLORIDyn.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ufechner7/FLORIDyn.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Introduction
A dynamic wind farm simulation software, translated from https://github.com/TUDelft-DataDrivenControl/FLORIDyn_Matlab, which was written by Marcus Becker.

His code uses one engineering model from the quasi static wind farm simulation software [FLORIS](https://github.com/NREL/floris), developed by NREL.

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

