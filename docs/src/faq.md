# Frequently asked questions FAQ

## Why should I use Julia?
**Short answer:** Speed (this code is 3-5 times faster than Matlab and 10-50 times faster than Python), which is very useful for example for solving optimization problems.
**Long answer:** Read [Why am I using Julia](https://ufechner7.github.io/2022/08/13/why-julia.html).

## How can I start Julia?
There are different options. Suggested way:
- Use a [Bash](https://www.w3schools.com/bash/index.php) terminal. On Linux this is the default, on Windows Bash is included in [Git for Windows](https://gitforwindows.org/). You can launch Bash from VSCode from the menu with `View->Terminal`. There is a small drop-down menu at the top left of the terminal window where you might have to select Bash if it is not the default.
- Basic method to launch Julia: Type `julia --project` in the Bash terminal.
- Improved method: Type `./bin/run_julia`. This script installs missing packages if needed and loads `Revise` and `FLORIDyn`.
- Expert method: Add the line `alias jl='bin/run_julia'` to your `.bashrc` file. Now you can start Julia by just typing `jl`.
Do NOT use the `run` button from VSCode to run Julia.

## Where can I find Julia packages?
If you need extra packages to solve your tasks, look at: [https://juliahub.com/ui/Packages](https://juliahub.com/ui/Packages) . They have a good search function. Don't ask AI, they often suggest outdated packages. If there are multiple packages for your problem, you can also ask on [Discourse](https://discourse.julialang.org/) for suggestions.

## What is a Julia environment?
A Julia environment is a folder that contains a `Project.toml` and a `Manifest.toml` file.
The `Project.toml` contains a list of the packages that are needed to run your code. It can also contain a list of the compatible package versions and more. It is important to use separate projects for each of your pieces of software that you develop, if you are not doing that (and use only the global environment), the things will break after some time. Further reading: [Working with Julia projects](https://ufechner7.github.io/2022/08/16/julia-projects.html).

### What should be in the global environment?
Start Julia with `julia`. On my PC I have:
```
julia> using Pkg; Pkg.status()
Status `~/.julia/environments/v1.11/Project.toml`
  [23c2ee80] ControlPlots v0.2.7
  [5903a43b] Infiltrator v1.9.2
  [16fef848] LiveServer v1.5.0
  [295af30f] Revise v3.8.1
  [0c614874] TerminalPager v0.6.4
  [1e6cf692] TestEnv v1.102.1
  [21f18d07] Timers v0.1.5
```
You can have a few more packages in there. But if you have 20 packages in your global environment you did something wrong.

### What should be in my project environment?
If you followed the developer guide and launched Julia with `jl` or `julia --project`, you should see something like:
```
julia> using Pkg; Pkg.status()
Project FLORIDyn v0.1.0
Status `~/repos/FLORIDyn.jl/Project.toml`
  [336ed68f] CSV v0.10.15
  [a93c6f00] DataFrames v1.7.0
  [8bb1440f] DelimitedFiles v1.9.1
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
  [9a3f8284] Random v1.11.0
Info Packages marked with ⌅ have new versions available but compatibility constraints restrict them from upgrading. To see why use `status --outdated`
``` 

## Which operating systems are supported?
Linux, Windows and Mac are supported. In some of the examples you might have to replace `/` with `\\` on Windows. If you still have a problem, create an issue on [https://github.com/ufechner7/FLORIDyn.jl/issues](https://github.com/ufechner7/FLORIDyn.jl/issues).

## How to create a new simulation?
- in the folder data, copy `2021_9T_Data.yaml` to a new file in the same folder. 
- copy the subfolder `2021_9T_Data` with the `.csv` files to a new folder with the same name as your new `.yaml` file (but without the suffix `.yaml`)
- copy the file `examples/main.jl` to a new file
- edit the line `settings_file = "data/2021_9T_Data.yaml"` to match the name of your new configuration file
- edit the new input files to match your test case
- run your new example using "include("examples/<new_main.jl>")"

## Can I use other plotting packages?
[ControlPlots.jl](https://github.com/aenarete/ControlPlots.jl) is an easy-to-use, powerful plotting package, suitable for teaching. Furthermore it is based on Matplotlib, so you can leverage your knowledge of Matplotlib if you have used Python before. [ControlPlots.jl](https://github.com/aenarete/ControlPlots.jl) exports the variable `plt`, and you can use any Matplotlib command by using the prefix `plt`.

You can use other plotting packages, but then you have to adapt the plotting scripts that can be found [here](https://github.com/ufechner7/FLORIDyn.jl/blob/main/src/visualisation/plot_flowfield.jl) yourself. Pull requests to support other plotting packages like `Plots.jl` or `Makie.jl` are welcome.

## Can I use this package with Python?
You can easily use Python packages in your own Julia project. I would suggest [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) to do that. You can use [Conda.jl](https://github.com/JuliaPy/Conda.jl) to install the required Python packages. They become part of your Julia environment and can be managed by the Julia package manager.

The other way, to use `FLORIDyn.jl` from Python is - in theory - possible using [PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl). This was not yet tested. 

## I have a Problem. Where can I get help?
Ask your question at [Discourse](https://discourse.julialang.org/). Most of the times you will get an answer in 15 min. Half of the people who answer are scientists.