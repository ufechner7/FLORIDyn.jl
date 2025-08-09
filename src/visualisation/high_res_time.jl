# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# High-resolution timestamp functions for FLORIDyn.jl
# Functions to get current time with microsecond and nanosecond resolution

"""
    now_microseconds()::String

Returns current timestamp as a string with microsecond resolution.
Format: "YYYY-mm-ddTHH-MM-SS.uuuuuu"

# Examples
```julia
julia> now_microseconds()
"2025-08-08T16-58-55.494911"
```
"""
function now_microseconds()::String
    # Get current time with nanosecond precision
    ns_time = time_ns()
    
    # Convert to DateTime (this gives us the base time)
    now_time = now()
    
    # Extract microseconds from nanoseconds
    # time_ns() gives nanoseconds since Unix epoch
    # We want the fractional second part in microseconds
    microseconds = div(mod(ns_time, 1_000_000_000), 1000)
    
    # Format the base datetime
    base_str = Dates.format(now_time, "yyyy-mm-ddTHH-MM-SS")
    
    # Add microseconds (6 digits with leading zeros)
    return string(base_str, ".", lpad(microseconds, 6, '0'))
end

"""
    now_nanoseconds()::String

Returns current timestamp as a string with nanosecond resolution.
Format: "YYYY-mm-ddTHH-MM-SS.nnnnnnnnn"

# Examples
```julia
julia> now_nanoseconds()
"2025-08-08T16-58-55.494911123"
```
"""
function now_nanoseconds()::String
    # Get current time with nanosecond precision
    ns_time = time_ns()
    
    # Convert to DateTime (this gives us the base time)
    now_time = now()
    
    # Extract nanoseconds from the fractional second part
    nanoseconds = mod(ns_time, 1_000_000_000)
    
    # Format the base datetime
    base_str = Dates.format(now_time, "yyyy-mm-ddTHH-MM-SS")

    # Add nanoseconds (9 digits with leading zeros)
    return string(base_str, ".", lpad(nanoseconds, 9, '0'))
end

"""
    precise_now()::String

Alias for now_microseconds() - returns current timestamp with microsecond resolution.
"""
precise_now() = now_microseconds()

"""
    unique_name()::String

Creates a unique directory name for storing simulation results.

This function generates a unique directory name by combining the prefix "floridyn_run_"
with a high-resolution timestamp. The timestamp includes microsecond precision to ensure
uniqueness even when multiple simulations are started in quick succession.

# Returns
- `String`: A unique directory name in the format `"floridyn_run_YYYY-mm-ddTHH-MM-SS.uuuuuu"`

# Examples
```julia
julia> unique_name()
"floridyn_run_2025-08-08T16-58-55.494911"

julia> unique_name()
"floridyn_run_2025-08-08T16-58-55.495123"
```

# Notes
- The generated name is suitable for use as a directory name on all operating systems
- Uses hyphens instead of colons in the time portion for filesystem compatibility
- Microsecond precision ensures uniqueness for rapid successive calls
- Used to create separate output directories for each run

# See Also
- [`now_microseconds`](@ref): The underlying timestamp function used
- [`Vis`](@ref): Visualization settings that may use unique names for output directories
"""
function unique_name()
    return "floridyn_run_" * now_microseconds()
end
