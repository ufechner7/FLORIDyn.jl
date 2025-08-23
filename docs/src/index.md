```@meta
CurrentModule = FLORIDyn
Description = "Dynamic wind farm simulation software"
```

# FLORIDyn

## Introduction
A dynamic wind farm simulation software, translated from [FLORIDyn_Matlab](https://github.com/TUDelft-DataDrivenControl/FLORIDyn_Matlab), which was written by Marcus Becker.
The code uses the Gaussian wake model derived in [3].

## Model features
- Simulate wind farms dynamically at a low computational cost
- Estimate the power generated, added turbulence, and wake-induced losses.
- Apply heterogeneous and time-varying wind speeds and directions
- Test different modeling approaches

![Flow Field](flowfield.png)

## Status:
Many examples work, some selectable via a menu:
```julia
include("examples/menu.jl")
```
The other examples can be executed directly using the `include` statement. Often, more than 30x the performance of
the Matlab version can be achieved.

A Python version of FLORIDyn is available at [https://github.com/TUDelft-DataDrivenControl/OFF](https://github.com/TUDelft-DataDrivenControl/OFF) .

## Installation
Install [Julia 1.11](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html), if you haven't already. Julia 1.10 is still supported. On Linux, make sure that Python3 and Matplotlib are installed:
```
sudo apt install python3-matplotlib
```
 
Make sure that `ControlPlots.jl` works as explained [here](https://github.com/aenarete/ControlPlots.jl?tab=readme-ov-file#installation).


Before installing this software it is suggested to create a new project, for example like this:
```bash
mkdir test
cd test
julia --project=.
```
Don't forget to type the `dot` at the end of the last command.

Then add FLORIDyn from  Julia's package manager, by typing:
```julia
using Pkg
pkg"add https://github.com/ufechner7/FLORIDyn.jl"
``` 
at the Julia prompt. You can run the unit tests with the command:
```julia
pkg"test FLORIDyn"
```
You can install the examples using the following command:
```julia
using FLORIDyn
install_examples()
```
If you now quit Julia with <ctrl><d> and restart it with
```bash
./bin/run_julia
```
then you can get the example menu by typing:
```julia
menu()
```
You can select any of the examples with the `<UP>` and `<DOWN>` keys, and then press `<ENTER>`.

## Installation using GIT
For developers, follow the [developer notes](https://ufechner7.github.io/FLORIDyn.jl/dev/developer/).

## References
Citation of the FLORIDyn model:\
[1] FLORIDyn - A dynamic and flexible framework for real-time wind farm control, M. Becker, D. Allaerts, J.W. van Wingerden, 2022, [doi:10.1088/1742-6596/2265/3/032103](http://doi.org/10.1088/1742-6596/2265/3/032103)

Used FLORIS model:\
[2] Experimental and theoretical study of wind turbine wakes in yawed conditions, M. Bastankhah, F. Porté-Agel, 2020, [doi:10.1017/jfm.2016.595](http://doi.org/10.1017/jfm.2016.595)

Gaussian wake model:\
[3] Experimental and theoretical study of wind turbine wakes in yawed conditions", M. Bastankhah, F. Porté-Agel, 2016, Journal of Fluid Mechanics. 2016;806:506-541. [doi:10.1017/jfm.2016.595](http://doi.org/10.1017/jfm.2016.595)

Additional references for smaller subcomponents can be found in the code or in the related publications.
