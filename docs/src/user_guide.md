# User Guide

This page contains information for users of FLORIDyn.jl.

## Installation of FLORIDyn
Installation is explained in the section [Installation](@ref).

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
Press `]` and then type `st` to see the active environment and its packages:
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
Remove unused packages with `rm <PACKAGE_NAME>`. Update outdated packages with `up`. Leave package mode with Backspace (←) and quit Julia with `Ctrl+D`.

### Launching Julia in the local environment
#### Standard way
Open a terminal in the project folder and type `julia --project`. Enter package mode (`]`) and you should see:
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

#### Advanced way
Starting Julia via a helper script can be faster and more consistent. Two scripts are provided:
- `./bin/run_julia`: Checks for missing dependencies and starts Julia with 11 threads (adjust to your CPU’s fast cores).
- `./bin/run_julia2`: Same check, starts Julia single-threaded.

Create aliases:

Linux / Windows (bash):
```bash
alias jl='./bin/run_julia'
alias jl2='./bin/run_julia2'
source ~/.bashrc
```

macOS (zsh or bash):
```bash
alias jl='./bin/run_julia'
alias jl2='./bin/run_julia2'
source ~/.zshrc          # for zsh
# or
source ~/.bash_profile   # for bash
```

Now run multi-threaded with `jl` or single-threaded with `jl2`.

## Multithreading
Modern CPUs have multiple cores. Threading lets independent loop iterations run in parallel using `Threads.@threads`. Only use it when iterations don’t depend on previous results. Threads share memory; immutable / read-only inputs can be shared; each thread needs separate mutable buffers. Writing to disjoint slices of an output array is thread-safe without locks.

If parallel functions allocate heavily, garbage collection pressure grows with thread count, limiting scaling unless allocation is reduced.

Currently only [getMeasurements](@ref) uses multithreading; it can speed up the simulation (with flow-field calculation) by about 4–5×.

## Multitasking
Running a simulation with online visualization in the same process is problematic: GUI updates cost time and the visualization library isn’t thread-safe. Therefore a second process handles visualization. It is started up-front (adds ~5 seconds time-to-first-plot). The line `@spawnat 2 Main.rmt_plot_flow_field(wf, X, Y, Z, vis; msr=msr)` runs the remote plotting function on process 2 and transfers the needed data efficiently.

## Running the examples
Start Julia with one of the scripts, then call `menu()`:
```
julia> menu()

Choose function to execute or `q` to quit:
 > select_project();                  print(CLEAR_SCR)
   select_measurement();              print(CLEAR_SCR)
   "plot_flow_field";                 PLT=1; include("main.jl")
   "plot_measurements";               PLT=4; include("main.jl")
   "plot_measurements_lineplot";      PLT=5; include("main.jl")
   "flow_field_vel_reduction_online"; PLT=6; include("main.jl")
   "create_video_from_saved_frames";  PLT=7; include("main.jl")
   "run_all_visualisations";          include("main_all.jl")
   "read_results";                    include("read_results.jl")
   "plot_wind_direction";             include("plot_wind_dir.jl")
   "play_videos";                     include("play_video.jl")
   open_documentation()
   quit
```
First select a project and a measurement. Then pick a visualization with the cursor keys and press Enter. Some examples aren’t in the menu; run them with `include("examples/<EXAMPLE>.jl")`. 
## Creating a new project
`data/projects.yaml`:
```yaml
projects:
  - project:
      name: 2021_9T_Data
      description: A reference simulation with 9 turbines
      settings: 2021_9T_Data.yaml
      vis: vis_default.yaml
  - project:
      name: 2021_9T_Data_full
      description: 9 turbines, no visualization, full output incl. 3 videos
      settings: 2021_9T_Data.yaml
      vis: vis_full_9T.yaml
  - project:
      name: 2021_54T_NordseeOne
      description: A reference simulation with 54 turbines
      settings: 2021_54T_NordseeOne.yaml
      vis: vis_54T.yaml
```
Add a new project by combining a `settings` file and a `vis` (visualization) file. To create a custom visualization file, copy an existing `vis_*.yaml`, rename it, adjust its content, and add a project entry. Then select it via `menu()`.

## Running a custom simulation
Steps:
1. Copy an existing settings YAML file in `data/` and give it a meaningful name.  
   The file `turbine_specs.yaml` need not be copied; just append new turbine definitions if needed.
2. Modify the settings YAML as required (follow inline comments).
3. Create a subfolder (named like the settings file, without extension) and copy all required CSV files.
4. Generate or update the CSV files.
5. Add a new project entry (see previous section).
6. Use `menu()` to run the simulation and visualize results.

### Optional
1. Copy `main_mini.jl` and adapt it.
2. Run it with `include("examples/<my_script>.jl")`.