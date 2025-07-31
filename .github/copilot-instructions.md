# Copilot instructions for Julia development

## Project context
This project is developed in Julia for scientific computing.

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

