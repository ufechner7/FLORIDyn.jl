# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Contributing to FLORIDyn.jl

Thank you for your interest in contributing to FLORIDyn.jl! This document provides guidelines for contributing to the project.

## Table of Contents

- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Contribution Guidelines](#contribution-guidelines)
- [Code Style](#code-style)
- [Testing](#testing)
- [Documentation](#documentation)
- [Pull Request Process](#pull-request-process)
- [Issue Reporting](#issue-reporting)
- [Performance Guidelines](#performance-guidelines)
- [License](#license)

## Getting Started

FLORIDyn.jl is a dynamic wind farm simulation software written in Julia. Before contributing, please:

1. Read the [documentation](https://ufechner7.github.io/FLORIDyn.jl/dev/)
2. Check existing [issues](https://github.com/ufechner7/FLORIDyn.jl/issues) and [pull requests](https://github.com/ufechner7/FLORIDyn.jl/pulls)
3. Review the [developer notes](https://ufechner7.github.io/FLORIDyn.jl/dev/developer/)

## Development Setup

### Prerequisites

1. **Julia 1.10 or higher** (Julia 1.11 recommended)
   - Install via [JuliaUp](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html)

2. **System dependencies** (Linux):
   ```bash
   sudo apt install python3-matplotlib
   ```

3. **ControlPlots.jl setup**
   - Follow the installation guide [here](https://github.com/aenarete/ControlPlots.jl?tab=readme-ov-file#installation)

### Clone and Setup

1. Fork the repository on GitHub
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/FLORIDyn.jl.git
   cd FLORIDyn.jl
   ```

3. Set up the development environment:
   ```bash
   ./bin/update
   ```

4. Set up convenient aliases (add to your `.bashrc`):
   ```bash
   alias jl='./bin/run_julia'
   alias jl2='./bin/run_julia2'
   ```

### Enable Revise for Development

For faster development iteration:
```bash
source ./bin/revise_on  # Enable Revise
source ./bin/revise_off # Disable Revise
```

## Contribution Guidelines

### Types of Contributions

We welcome various types of contributions:

- **Bug fixes**: Fix identified issues
- **New features**: Implement new wind farm modeling capabilities
- **Performance improvements**: Optimize computational efficiency
- **Documentation**: Improve docs, examples, or tutorials
- **Tests**: Add or improve test coverage
- **Examples**: Create new demonstration scripts

### Before You Start

1. **Discuss major changes**: Open an issue to discuss significant features or changes
2. **Check for duplicates**: Ensure your contribution doesn't duplicate existing work
3. **Start small**: Begin with small contributions to familiarize yourself with the codebase

## Code Style

### Julia Style Guidelines

Follow these coding standards:

- **Indentation**: Use 4 spaces (no tabs)
- **Function style**: Write reusable functions with explicit parameters
- **Conditionals**: Avoid parentheses around conditions
  ```julia
  # Good
  if x > 0
      return x
  end
  
  # Avoid
  if (x > 0)
      return x
  end
  ```
- **Type annotations**: Avoid overly specific types unless necessary
- **Broadcasting**: Use `.` for element-wise operations
- **Multiple dispatch**: Leverage Julia's multiple dispatch system

### Performance Guidelines

- **Functions over scripts**: Place performance-critical code inside functions
- **Avoid global variables**: Minimize untyped global variables
- **Memory efficiency**: Reuse objects when possible to minimize allocations
- **Pre-allocation**: Use pre-allocated buffers for performance-critical loops

### Documentation Standards

- **Docstrings**: Include docstrings for all public functions and modules
- **Type links**: Add links to type definitions in docstrings for parameters with type assertions
- **Examples**: Provide usage examples in docstrings
- **Conciseness**: Be concise but comprehensive

Example docstring format:
```julia
"""
    myFunction(param1::Type1, param2::Type2) -> ReturnType

Brief description of what the function does.

# Arguments
- `param1::Type1`: Description of param1. See: [`Type1`](@ref)
- `param2::Type2`: Description of param2

# Returns
- Description of return value

# Examples
```julia
result = myFunction(value1, value2)
```
"""
```

## Testing

### Running Tests

Test your changes thoroughly:

```julia
# Run all tests
using Pkg
Pkg.test()

# Run specific test files
include("test/test_floris.jl")

# Interactive test menu
include("test/runtest.jl")
```

### Test Requirements

- **Single and multi-threaded**: Test both execution modes:
  ```bash
  jl   # Single-threaded
  ]test
  exit()
  jl2  # Multi-threaded
  ]test
  exit()
  ```

- **New functionality**: Add tests for any new features or bug fixes
- **Edge cases**: Test boundary conditions and error cases
- **Performance**: Ensure new code doesn't significantly degrade performance

### Test Structure

- Place tests in the `test/` directory
- Follow the naming convention: `test_<module>.jl`
- Use Julia's `@testset` macros for organization
- Include both unit tests and integration tests

## Documentation

### Building Documentation

Generate documentation locally:
```julia
include("scripts/build_docu.jl")
```

### Documentation Requirements

- **New features**: Document all new public functions and types
- **Examples**: Provide working examples for new functionality
- **Developer notes**: Update developer documentation for significant changes
- **Changelog**: Mention noteworthy changes (maintainers will handle formal changelog)

## Pull Request Process

### Before Submitting

1. **Test thoroughly**: Ensure all tests pass
2. **Check performance**: Verify no significant performance regressions
3. **Update documentation**: Add/update relevant documentation
4. **Code style**: Follow the established coding conventions
5. **Commit messages**: Write clear, descriptive commit messages

### PR Submission

1. **Create a branch**: Use a descriptive branch name
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make commits**: Make logical, well-documented commits
   ```bash
   git add .
   git commit -m "Add feature: brief description of changes"
   ```

3. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

4. **Open a Pull Request**: Include:
   - Clear description of changes
   - Motivation for the changes
   - Testing performed
   - Any breaking changes

### PR Review Process

- Maintainers will review your PR
- Address any requested changes
- Once approved, maintainers will merge the PR

## Issue Reporting

### Bug Reports

When reporting bugs, include:

- **Julia version**: Output of `versioninfo()`
- **FLORIDyn.jl version**: Package version being used
- **System information**: OS, architecture
- **Minimal example**: Code that reproduces the issue
- **Expected vs actual behavior**: Clear description of the problem
- **Stack trace**: If applicable, include the full error message

### Feature Requests

For feature requests, provide:

- **Use case**: Why is this feature needed?
- **Proposed solution**: How should it work?
- **Alternatives considered**: What other approaches were considered?
- **Implementation ideas**: Any thoughts on implementation

## Performance Guidelines

### Memory Allocation

Keep allocations low for good performance:

- **GC time < 10%**: Garbage collection time should be less than 10% of total execution
- **Profile allocations**: Use the provided tools:
  ```bash
  ./bin/test_alloc    # Run allocation profiling
  ./bin/analyze_alloc # Analyze results
  ```

### Optimization Tips

- **Inner loops**: Focus optimization on innermost functions/loops
- **Type stability**: Ensure functions are type-stable
- **Pre-allocation**: Use pre-allocated buffers for repeated operations
- **Views**: Use `@views` to avoid unnecessary allocations

## Debugging

### Recommended Tools

Use [Infiltrator.jl](https://github.com/JuliaDebug/Infiltrator.jl) for debugging:

1. Launch debug session: `jl2`
2. Load package: `using FLORIDyn`
3. Add breakpoint: `Main.@infiltrate`
4. Run your code
5. Inspect variables at the `infil>` prompt

## License

By contributing to FLORIDyn.jl, you agree that your contributions will be licensed under the BSD-3-Clause license. See the [LICENSE](LICENSE) file for details.

## Questions?

- **Documentation**: Check the [developer notes](https://ufechner7.github.io/FLORIDyn.jl/dev/developer/)
- **Issues**: Open an issue for bugs or feature requests
- **Discussions**: Use [Discourse](https://discourse.julialang.org/) for general questions and discussions

Thank you for contributing to FLORIDyn.jl! 🌬️
