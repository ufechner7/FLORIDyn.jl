# User guide

This page contains information for users of FLORIDyn.jl.

## Installation
Currently, installation using git is recommended. Later this year installation as Julia package will be better tested. This requires the issue [Create a package extension for plotting](https://github.com/ufechner7/FLORIDyn.jl/issues/39)
to be resolved.

## Launching Julia
### Launching Julia in the global environment
Just type `julia`. 
```julia
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.11.6 (2025-07-09)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> 
```
If you now type `]` and then `st` you can see the name of the environment and which packages it contains:
```
(@v1.11) pkg> st
Status `~/.julia/environments/v1.11/Project.toml`
  [23c2ee80] ControlPlots v0.2.7
  [fab6aee4] DistributedNext v1.1.0
  [5903a43b] Infiltrator v1.9.2
  [16fef848] LiveServer v1.5.0
  [5f6e1e16] LocalCoverage v0.8.3
  [295af30f] Revise v3.8.1
  [0c614874] TerminalPager v0.6.4
  [1e6cf692] TestEnv v1.102.1
  [21f18d07] Timers v0.1.5
```
If there are unused packages in your global environment, then remove them by typing `rm <PACKAGE_NAME>`, if packages are outdated update them with `up`. You can leave the package manager mode by typing `<BACK>` and quit Julia by typing `<CTRL>+<D>`.

### Launching Julia in the local environment

## Multithreading

## Multitasking

## Running the examples

