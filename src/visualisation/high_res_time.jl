# High-resolution timestamp functions for FLORIDyn.jl
# Functions to get current time with microsecond and nanosecond resolution

# Note: Dates is a standard library, so it should be available

"""
    now_microseconds()::String

Returns current timestamp as a string with microsecond resolution.
Format: "YYYY-mm-ddTHH:MM:SS.μμμμμμ"

# Examples
```julia
julia> now_microseconds()
"2025-08-08T16:58:55.494911"
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
    microseconds = div(ns_time % 1_000_000_000, 1000)
    
    # Format the base datetime
    base_str = Dates.format(now_time, "yyyy-mm-ddTHH:MM:SS")
    
    # Add microseconds (6 digits with leading zeros)
    return string(base_str, ".", lpad(microseconds, 6, '0'))
end

"""
    now_nanoseconds()::String

Returns current timestamp as a string with nanosecond resolution.
Format: "YYYY-mm-ddTHH:MM:SS.nnnnnnnnn"

# Examples
```julia
julia> now_nanoseconds()
"2025-08-08T16:58:55.494911123"
```
"""
function now_nanoseconds()::String
    # Get current time with nanosecond precision
    ns_time = time_ns()
    
    # Convert to DateTime (this gives us the base time)
    now_time = now()
    
    # Extract nanoseconds from the fractional second part
    nanoseconds = ns_time % 1_000_000_000
    
    # Format the base datetime
    base_str = Dates.format(now_time, "yyyy-mm-ddTHH:MM:SS")
    
    # Add nanoseconds (9 digits with leading zeros)
    return string(base_str, ".", lpad(nanoseconds, 9, '0'))
end

"""
    precise_now()::String

Alias for now_microseconds() - returns current timestamp with microsecond resolution.
"""
precise_now() = now_microseconds()

# Export the functions if this file is included in a module
# export now_microseconds, now_nanoseconds, precise_now
