# FLORIDyn
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ufechner7.github.io/FLORIDyn.jl/dev)
[![DOI](https://zenodo.org/badge/DOI/17121846.svg)](https://doi.org/17121846)
[![Build Status](https://github.com/ufechner7/FLORIDyn.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ufechner7/FLORIDyn.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/ufechner7/FLORIDyn.jl/graph/badge.svg?token=O7wXT62VSR)](https://codecov.io/gh/ufechner7/FLORIDyn.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Introduction
A dynamic wind farm simulation software, translated from [FLORIDyn_Matlab](https://github.com/TUDelft-DataDrivenControl/FLORIDyn_Matlab) [\[1\]](#ref1), which was written by Marcus Becker.
The code uses the Gaussian wake model derived in [\[3\]](#ref3).

## Model features
- Simulate wind farms dynamically at a low computational cost
- Estimate the power generated, added turbulence, and wake-induced losses.
- Apply heterogeneous and time-varying wind speeds and directions
- Test different modeling approaches

<p align="center"><img src="https://github.com/ufechner7/FLORIDyn.jl/blob/main/docs/src/flowfield.png?raw=true" width="500" /></p>

## Status
All examples work, most key examples are selectable via a menu:
```julia
include("examples/menu.jl")
```
The other examples can be executed directly using the `include` statement. Often, more than 30x the performance of the Matlab version can be achieved. Currently, only `IterateOPs_basic` is implemented.

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
pkg"add FLORIDyn"
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
You can select any of the examples with the \<UP\> and \<DOWN\> keys, and then press \<ENTER\>.

## Installation using GIT
For developers, follow the [developer notes](https://ufechner7.github.io/FLORIDyn.jl/dev/developer/).

## Documentation
The documentation is available [here](https://ufechner7.github.io/FLORIDyn.jl/dev/).

## License
This project is licensed under the  `BSD-3-Clause`. The documentation is licensed under the `CC-BY-4.0 License`. Please see the below `Copyright notice` in association with the licenses that can be found in the file [LICENSE](LICENSE) in this folder.

## Copyright notice
Technische Universiteit Delft hereby disclaims all copyright interest in the package “FLORIDyn.jl” (dynamic wind farm simulation) written by the Author(s).

Prof.dr. H.G.C. (Henri) Werij, Dean of Aerospace Engineering, Technische Universiteit Delft.

See the copyright notices in the source files, and the list of authors in [AUTHORS.md](AUTHORS.md).

## Sponsors
This research has been partly funded by Rijksdienst voor Ondernemend (RVO) Nederland under contract HEP24-03681024 through the EuroWindWakes project which is a Horizon Europe Partnership.

## References
<a id="ref1"></a>
Citation of the FLORIDyn model:  
[1] FLORIDyn - A dynamic and flexible framework for real-time wind farm control, M. Becker, D. Allaerts, J.W. van Wingerden, 2022, [DOI 10.1088/1742-6596/2265/3/032103](http://doi.org/10.1088/1742-6596/2265/3/032103)

Used FLORIS model:  
[2] Experimental and theoretical study of wind turbine wakes in yawed conditions, M. Bastankhah, F. Porté-Agel, 2020, [DOI 10.1017/jfm.2016.595](http://doi.org/10.1017/jfm.2016.595)

<a id="ref3"></a>
Gaussian wake model:  
[3] Experimental and theoretical study of wind turbine wakes in yawed conditions, M. Bastankhah, F. Porté-Agel, 2016, Journal of Fluid Mechanics 806:506-541. [DOI 10.1017/jfm.2016.595](http://doi.org/10.1017/jfm.2016.595)

Additional references for smaller subcomponents can be found in the code or in the related publications.

## Publications that use FLORIDyn
- Ensemble-Based Flow Field Estimation Using the Dynamic Wind Farm Model FLORIDyn, M. Becker, D. Allaerts, J.W. van Wingerden, 2022, [DOI 10.3390/en15228589](http://doi.org/10.3390/en15228589)
- Wind pattern clustering of high frequent field measurements for dynamic wind farm flow control, M. Becker, D. Allaerts, J.W. van Wingerden, 2024, http://doi.org/10.1088/1742-6596/2767/3/032028 
- Sensitivity analysis and Bayesian calibration of a dynamic wind farm control model: FLORIDyn, V.V. Dighe, M. Becker, wf. Göçmen, B. Sanderse, J.W. van Wingerden, 2022, http://doi.org/10.1088/1742-6596/2265/2/022062
- Time-shifted cost function design for more efficient dynamic wind farm flow control, M. Becker, D. Allaerts, J.W. van Wingerden, 2024, http://doi.org/10.1109/CCTA60707.2024.10666535
- Suitability of Dynamic Wake Models for AEP Estimation: A Wind Farm-Scale Validation Study, M. Van der Straeten, http://resolver.tudelft.nl/uuid:f35617a2-2409-439b-8bc2-6334b807ce1f 
- Scaling DMD modes for modeling Dynamic Induction Control wakes in various wind speeds, J. Gutknecht, M. Becker, C. Muscari, wf. Lutz, J.W. van Wingerden, 2023, http://doi.org/10.1109/CCTA54093.2023.10252400
- Model predictive control of wakes for wind farm power tracking, A. Sterle, C.A. Hans, J. Raisch, 2024, http://doi.org/10.1088/1742-6596/2767/3/032005

