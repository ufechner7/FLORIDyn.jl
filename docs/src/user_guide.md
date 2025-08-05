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
Modern CPUs can have many cores. Threading is one way to make to split calculations such that the are calculated
in parallel on multiple cores. In Julia, the `@threads` macro can be used in front of a for loop to do calculations
in parallel. This only works if the iterations are independent of each other, it will not work if iteration `n+1`
depends on results calculated by iteration `n`. Threads - in contrast to processes - share the same memory. If the
input data is immutable and accessed in a read only way, the threads can use the same input data. Each thread needs
his own buffer for mutable data. If each thread writes to a unique, non-overlapping portion of an output array, 
this operation is thread-safe and does not require locks or other synchronization mechanisms.

If the functions that are executed in parallel allocate memory, then the pressure on the garbage collector increases
with the number of threads. This creates a practical limit for the number of threads that can be used in one process,
unless the parallel code is allocation free.

In this simulation software, currently only the function [getMeasurementsP](@ref) uses multithreading. It can increase 
the performance of the simulation with flow field calculation by a factor of four to five.

## Multitasking
Running a simulation with online visualization is problematic, because firstly updating the GUI costs time, 
and secondly the currently used visualization library is not thread safe. Therefore we use a second process for
visualization. The second process is started at the beginning, which increases the time-to-first-plot by about
five seconds. The line ` @spawnat 2 Main.rmt_plot_flow_field(wf, X, Y, Z, vis; msr=msr)` calls the remote plotting
function on the second process and transfers the required data. This is safe and very fast.

## Running the examples

To run the examples, launch Julia with one of the start scripts and then type `menu()`. You should see:
```
julia> menu()

Choose function to execute or `q` to quit: 
 > flow_field_vel_reduction        = PLT=1; include("main.jl")
   flow_field_added_turbulence     = PLT=2; include("main.jl")
   flow_field_eff_wind_speed       = PLT=3; include("main.jl")
   plot_measurements_              = PLT=4; include("main.jl")
   plot_measurements_lineplot      = PLT=5; include("main.jl")
   flow_field_vel_reduction_online = PLT=6; include("main.jl")
   create_video_from_saved_frames  = PLT=7; include("main.jl")
   play_videos                     = include("play_video.jl")
   open_documentation()
   quit
```
or similar. You can execute any of the examples by selecting one of them with the cursor keys and then press 
<ENTER>. There might be additional examples that are not yet integrated in the menu. You can execute them
with `include("examples/<MY_EXAMPLE.jl>")`. Some examples require that Julia runs in single-threaded mode.
If you want to run such an example, start Julia with `jl2`.

## Running a custom simulation
To run your own simulation you need to follow these steps:
1. Create a copy of an existing YAML file in the data folder. Give it a good name.
2. Modify the custom YAML file according to your needs, following the comments in the YAML file.
3. Create a subfolder with the name of the custom YAML file and copy all required CSV files.
4. Update/ generate the CSV files according to your needs.
5. Copy the script `main_mini.jl` and adapt it according to your needs.
6. Run your new script using the command `include("examples/<my_script.jl>")`