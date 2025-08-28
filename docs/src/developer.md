# Developer notes

This page contains information for developers who want to contribute to or extend FLORIDyn.jl.

## Project structure

The FLORIDyn.jl package is organized into several modules:

- `src/` - Main source code
  - `controller/`    - Controller implementations
  - `correction/`    - Wake correction models
  - `floridyn_cl/`   - Main simulation loop
  - `floris/`        - FLORIS model implementations
  - `windfield/`     - Wind field modeling
  - `visualisation/` - Plotting and helper functions
- `test/` - Test suite
- `examples/` - Example scripts
- `examples_dev/` - Examples for developers
- `docs/`  - Documentation source
- `data/`  - Example data files and configuration settings
- `video/` - Intermediate PNG output files
- `out/`   - Final PNG and MP4 output files, one folder per test run
- `bin/`   - Bash scripts to start Julia, analyze allocation and statistics

## Development workflow

### Setting up the development environment

1. On Linux, make sure that Python3 and Matplotlib are installed:
```
sudo apt install python3-matplotlib
```
 
Make sure that `ControlPlots.jl` works as explained [here](https://github.com/aenarete/ControlPlots.jl?tab=readme-ov-file#installation).

2. Clone the repository:
   ```bash
   git clone https://github.com/ufechner7/FLORIDyn.jl.git
   cd FLORIDyn.jl
   ```

3. Activate the project environment and update the packages:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.up()
   ```

### Running tests
To run the full test suite:
```julia
using Pkg
Pkg.test()
```

To run specific tests:
```julia
include("test/test_floris.jl")
```

You get a menu with the available test sets and benchmarks using:
```julia
include("test/runtest.jl")
```
Be aware that the tests should be run both single threaded and multithreaded. So for a full test do:
```bash
jl
]test
exit()
jl2
]test
exit()
```
**Make sure you have followed the steps in the section [Advanced way](@ref) of the user guide.**

### Using Revise
Using [Revise](https://github.com/timholy/Revise.jl) allows you to edit and run code that is part of `FLORIDyn.jl` without restarting Julia.
Limitation before Julia 1.12: If you modify structs, you still have to restart Julia for the changes to become into effect.

To enable `Revise` type
```bash
source ./bin/revise_on
```
to disable it
```bash
source ./bin/revise_off
```
Be aware that the environment variable `USE_REVISE` has only an impact if Julia is launched using `jl` or `jl2`.

### Debugging
Suggestion: Use [Infiltrator.jl](https://github.com/JuliaDebug/Infiltrator.jl) for debugging.

Add two lines to your `.bashrc` script (create one if it does not exist yet):
```bash
alias jl='./bin/run_julia'
alias jl2='./bin/run_julia2'
```
Install the packages `Infiltrator.jl` and `Revise.jl` in your global environment:
```bash
julia -e 'using Pkg; Pkg.add("Revise"); Pkg.add("Infiltrator")'
```
Now you can launch Julia by typing `jl`, and for debugging type `jl2`.

Debugging session:
1. bash> `jl2`
2. julia> `using FLORIDyn`
3. Add the line `Main.@infiltrate` at the location where you want to set a break point.
4. julia> `include("examples/main.jl")`
Now the program should run into your breakpoint. You should see the prompt:
```
infil> 
```
and by typing the name of any local or global variable you can inspect the content.
You can also execute any statement that fails and modify it until it works.

**Important:** When executing STEP 2, the line `Main.@infiltrate` must not exist
or must be commented.

### Building documentation

To build the documentation locally:
```julia
include("scripts/build_docu.jl")
```
You can get an overview over the exported methods by running the script:
```julia
include("scripts/stats.jl")
```

## Checking the memory allocations
For a good performance, in particular when using multi-threading, the amount of memory allocations should be low.
As a general rule, if you use a benchmark script and the time that the garbage collector (GC) needs is less
than 10%, that should be good enough.

For fixing allocations it is most of the time sufficient to fix the most inner functions or loops. To find them,
two Bash scripts are provided: `test_alloc` and `analyze_alloc` . The second script has an output like this:
```
ufechner@ufryzen:~/repos/FLORIDyn.jl/bin$ ./analyze_alloc 
Found 19 .mem files
Analyzing memory allocations...

=== TOP 20 MEMORY ALLOCATIONS ===

190.67 MB    ./windfield/windfield_shear.jl.61012:74 shear = z_norm .^ wind_shear.alpha
93.94 MB     ./floridyn_cl/floridyn_cl.jl.61012:492 runFLORIS!(
5.57 MB      ./floridyn_cl/floridyn_cl.jl.61012:604 T_addedTI = sqrt(sum(T_aTI_arr .^ 2))
4.48 MB      ./floridyn_cl/floridyn_cl.jl.61012:681 wf.Weight[iT] = wf.Weight[iT] ./ wS
3.48 MB      ./visualisation/calc_flowfield.jl.61012:246 GP.posNac[end, :] = [0.0, 0.0, zh] # Update grid point nacelle position
3.48 MB      ./visualisation/calc_flowfield.jl.61012:245 GP.posBase[end, :] = [xGP, yGP, 0.0] # Update grid point base position
269.10 KB    ./floridyn_cl/iterate.jl.61012:50 return IterateOPsBuffers(
168.77 KB    ./settings.jl.61012:654 diff = sum(
164.61 KB    ./visualisation/plot_flowfield.jl.61012:434 rm(joinpath("video", file))
128.70 KB    ./floris/gaussian.jl.61012:180 Ct = calcCt(states_t[:, 1], states_t[:, 2])
127.60 KB    ./floris/gaussian.jl.61012:170 RPs = Matrix{Float64}(undef, n, 3)
85.78 KB     ./floris/gaussian.jl.61012:181 yaw = -deg2rad.(states_t[:, 2])
84.46 KB     ./floris/gaussian.jl.61012:264 states_op = copy(wf.States_OP)
...

=== ALLOCATION SUMMARY BY FILE ===

190.67 MB    1    lines ./windfield/windfield_shear.jl.61012
103.99 MB    15   lines ./floridyn_cl/floridyn_cl.jl.61012
6.96 MB      3    lines ./visualisation/calc_flowfield.jl.61012
980.55 KB    32   lines ./floris/gaussian.jl.61012
287.91 KB    29   lines ./visualisation/plot_flowfield.jl.61012
269.10 KB    1    lines ./floridyn_cl/iterate.jl.61012
191.23 KB    25   lines ./floridyn_cl/prepare_simulation.jl.61012
190.89 KB    28   lines ./settings.jl.61012
20.30 KB     20   lines ./FLORIDyn.jl.61012
...
```
So you can clearly see where the largest allocations happen and refactor those functions or lines.

## Code style and conventions

### Naming conventions
- Use descriptive variable names
- Follow Julia naming conventions (lowercase with underscores for functions and variables)
- Type names should use CamelCase
- Constants should be ALL_CAPS

### Code organization
- Keep functions focused and small
- Use meaningful docstrings for all exported functions
- Include type annotations where helpful for clarity
- Follow the existing module structure

### Testing
- Write tests for new functionality
- Ensure all tests pass before submitting pull requests
- Include edge cases in test coverage
- Use descriptive test names

## Contributing

### Pull requests
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Update documentation if needed
7. Submit a pull request

### Issues
When reporting issues, please include:
- Julia version
- FLORIDyn.jl version
- Minimal working example
- Error messages and stack traces

## Architecture overview

### Core simulation loop
The main simulation is handled by the `floridyn_cl` module, which implements the time-stepping algorithm for wake evolution.

### Wake models
The package supports multiple wake models through abstract types:
- `VelModel` - Velocity deficit models
- `DirModel` - Wake deflection models
- `TurbulenceModel` - Turbulence intensity models

### Settings system
Configuration is handled through YAML files that are parsed into Julia structs. See the Settings documentation for details.

## Performance considerations
- The simulation uses in-place operations where possible to minimize allocations
- Key loops are optimized for performance
- Consider using `@profile` and `BenchmarkTools.jl` when optimizing code
- Read the [performance tips] (https://docs.julialang.org/en/v1/manual/performance-tips/)
