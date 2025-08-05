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
#### The standard way
Launch a terminal in the correct folder, and then type `julia --project`. If you now enter the package manager mode with `]`, you should see:
```
(FLORIDyn) pkg> st
Project FLORIDyn v0.1.0
Status `~/repos/FLORIDyn.jl/Project.toml`
  [336ed68f] CSV v0.10.15
  [a93c6f00] DataFrames v1.7.0
  [8bb1440f] DelimitedFiles v1.9.1
  [fab6aee4] DistributedNext v1.1.0
  [ffbed154] DocStringExtensions v0.9.5
  [a98d9a8b] Interpolations v0.16.1
  [033835bb] JLD2 v0.5.15
  [b964fa9f] LaTeXStrings v1.4.0
  [d96e819e] Parameters v0.12.3
⌅ [aea7be01] PrecompileTools v1.2.1
  [90137ffa] StaticArrays v1.9.14
  [10745b16] Statistics v1.11.1
  [7c3b921d] StructMapping v0.2.3
  [ddb6d928] YAML v0.4.14
  [37e2e46d] LinearAlgebra v1.11.0
  [56ddb016] Logging v1.11.0
  [44cfe95a] Pkg v1.11.0
  [de0858da] Printf v1.11.0
  [9a3f8284] Random v1.11.0
Info Packages marked with ⌅ have new versions available but compatibility constraints restrict them from upgrading. To see why use `status --outdated`
```

#### The advanced way
It is faster and gives better results if you start Julia using a script. Currently, two user scripts are provided:
- The script `./bin/run_julia` checks for missing dependencies and starts Julia with 11 threads. Adapt the script such that this number matches the number of fast cores of your CPU.
- The script `./bin/run_julia2` checks for missing dependencies and starts Julia with one thread.

To make these scripts easier to use, you can create bash aliases. 

**On Linux or Windows:**
Add the following lines to your `~/.bashrc` file:

```bash
alias jl='./bin/run_julia'
alias jl2='./bin/run_julia2'
```

After adding these aliases, reload your bash configuration:
```bash
source ~/.bashrc
```

**On Mac:**
Add the following lines to your `~/.zshrc` file (for zsh, which is the default shell on modern macOS) or `~/.bash_profile` file (for bash):

```bash
alias jl='./bin/run_julia'
alias jl2='./bin/run_julia2'
```

After adding these aliases, reload your shell configuration:
```bash
source ~/.zshrc    # for zsh
# or
source ~/.bash_profile    # for bash
```

Now you can launch Julia with multithreading by simply typing `jl` or with single threading by typing `jl2`.


## Multithreading

## Multitasking

## Running the examples

