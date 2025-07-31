# Copilot instructions for Julia development

## Project context
This project is developed in Julia for scientific computing.

## Folder Structure
- `/src`: Contains the source code of this package.
  - `/controller`: Contains implementations of various controllers.
  - `/correction`: Contains wake correction models.
  - `/floridyn_cl`: Contains the main simulation loop.
  - `/floris`: Contains FLORIS model implementations.
  - `/windfield`: Contains wind field modeling code.
  - `/visualisation`: Contains plotting and helper functions.
- `/test`: Contains the test suite for the package.
- `/examples`: Contains example scripts demonstrating package usage.
- `/docs`: Contains the documentation source files.
- `/docs/src`: Contains the documentation source files, including `types.md` and `developer.md`.
- `/data`: Contains example data files and configuration settings.
- `/video`: Contains PNG and MP4 output files.
- `/bin`: Contains bash scripts to start Julia and for statistics.

## Coding style
- Use 4 spaces for indentation.
- Write reusable functions instead of scripts; functions should take parameters explicitly.
- Avoid parentheses around conditions (e.g., write `if x > 0` rather than `if (x > 0)`).
- Avoid overly specific type annotations unless necessary.

## Performance
- Place performance-critical code inside functions.
- Avoid untyped global variables.
- Minimize memory allocations by reusing objects when possible.

## Idiomatic usage
- Use multiple dispatch to write flexible, generic functions.
- Use broadcasting (`.`) for element-wise operations.

## Documentation
- Include docstrings for functions and modules.
- Constructors must return an instance of their own type.
- Be concise
- For parameters that have a type assertion, add a link to the type definition to the docstring

## General
- Write clean, readable, and maintainable code.

