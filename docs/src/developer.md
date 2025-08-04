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
- `video/` - PNG and MP4 output files
- `bin/`   - Bash scripts to start Julia and for statistics

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
