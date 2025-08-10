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

"""
    delete_results(vis::Vis, n::Int=1; dry_run::Bool = false)

Delete the newest `n` directories starting with "floridyn_run" from both the visualization 
output directory and video directory.

This function provides comprehensive cleanup by removing the most recent floridyn_run directories
from both output and video paths simultaneously. It's particularly useful for removing failed runs, 
test runs, or managing disk space by keeping only the most relevant simulation outputs across
both directory types.

# Arguments
- `vis::Vis`: Visualization settings object containing both output and video path configurations
- `n::Int`: Number of newest directories to delete from each directory (must be positive, default: 1)  
- `dry_run::Bool`: Preview mode - shows what would be deleted without actually deleting (default: `false`)

# Returns
- `Vector{String}`: Absolute paths of directories that were deleted from both locations combined.
  Returns empty vector if no matching directories found or if `n ≤ 0`.

# Behavior
- **Directory Search**: Searches both `vis.output_path` and `vis.video_path` for directories matching `floridyn_run_*` pattern
- **Dual Cleanup**: Operates on both output and video directories in sequence
- **Sorting**: Sorts directories by modification time (newest first) within each directory
- **Selection**: Selects up to `n` newest directories for deletion from each location
- **Cleanup**: Automatically sets `vis.unique_folder = ""` before processing
- **Dry Run**: When `dry_run=true`, logs what would be deleted but performs no actual deletion

# Examples
```julia
# Create visualization settings
vis = Vis("data/vis_default.yaml")

# Delete the single newest floridyn_run directory from both output and video directories
deleted = delete_results(vis)
println("Deleted: ", length(deleted), " directories total")

# Delete the 3 newest floridyn_run directories from each location
deleted = delete_results(vis, 3)
println("Deleted directories: ", basename.(deleted))

# Preview what would be deleted from both directories
delete_results(vis, 5, dry_run=true)  # Shows info about 5 newest in each directory

# Check existing directories in both locations
vis.unique_folder = ""  
output_runs = find_floridyn_runs(vis.output_path)
video_runs = find_floridyn_runs(vis.video_path)
println("Found \$(length(output_runs)) in output, \$(length(video_runs)) in video")
```

# Error Handling
- **Invalid Input**: Warns and returns empty vector for non-positive `n`
- **Missing Directories**: Errors if either `vis.output_path` or `vis.video_path` doesn't exist
- **No Matches**: Info message and early return if no floridyn_run directories found in either location
- **Deletion Failures**: Individual directory deletion errors are logged but don't stop the process

# Important Notes
- **Dual Operation**: This function operates on BOTH output and video directories
- **Independent Processing**: Each directory is processed separately - `n` directories from output AND `n` directories from video
- **Early Return**: Function returns early if no directories are found in the first location checked
- **Path Creation**: Accessing `vis.output_path` and `vis.video_path` automatically creates these directories if they don't exist

# See Also
- [`Vis`](@ref): Visualization settings struct with output and video path configuration
- [`unique_name()`](@ref): Creates timestamped floridyn_run directories
- [`find_floridyn_runs()`](@ref): Lists existing floridyn_run directories in a given path
"""
function delete_results(vis::Vis, n::Int=1; dry_run::Bool = false)
    vis.unique_folder = ""
    deleted_dirs = String[]
    
    for directory in [vis.output_path, vis.video_path]
        if n <= 0
            @warn "Number of folders to delete must be positive. Got: $n"
            return String[]
        end
        
        if !isdir(directory)
            @error "Directory does not exist: $directory"
            return String[]
        end
        
        # Find all floridyn_run directories
        floridyn_dirs = find_floridyn_runs(directory)
        
        if isempty(floridyn_dirs)
            @info "No floridyn_run directories found in: $directory"
            return String[]
        end
        
        # Sort by modification time (newest first for deletion)
        sorted_dirs = sort(floridyn_dirs, by=d -> mtime(d), rev=true)
        
        # Determine how many to actually delete
        n_to_delete = min(n, length(sorted_dirs))
        dirs_to_delete = sorted_dirs[1:n_to_delete]
        
        if dry_run
            @info "DRY RUN - Would delete $n_to_delete newest floridyn_run directories:"
            for (i, dir) in enumerate(dirs_to_delete)
                dir_name = basename(dir)
                mod_time = Dates.unix2datetime(mtime(dir))
                @info "  $(i). $dir_name (modified: $mod_time)"
            end
            return dirs_to_delete
        end
        
        @info "Deleting $n_to_delete newest floridyn_run directories from: $directory"
        
        for (i, dir) in enumerate(dirs_to_delete)
            dir_name = basename(dir)
            try
                rm(dir, recursive=true)
                push!(deleted_dirs, dir)
                @info "  $(i)/$n_to_delete Deleted: $dir_name"
            catch e
                @error "Failed to delete $dir_name: $e"
            end
        end
        
        remaining_count = length(floridyn_dirs) - length(deleted_dirs)
        @info "Cleanup complete. Deleted $(length(deleted_dirs)) directories, $remaining_count remaining."
    end 
    return deleted_dirs
end

"""
    find_floridyn_runs(directory::String = pwd())

Find all directories starting with "floridyn_run" in the specified directory.

# Arguments  
- `directory::String`: Directory to search in (default: current working directory)

# Returns
- `Vector{String}`: List of full paths to floridyn_run directories, sorted by name

# Examples
```julia
# Find floridyn_run directories in current directory
dirs = find_floridyn_runs()

# Find in specific directory
dirs = find_floridyn_runs("out")

# Count how many exist
count = length(find_floridyn_runs())
```

# See Also
- [`delete_results`](@ref): Function that uses this to find directories to delete
- [`unique_name`](@ref): Function that creates these directories
"""
function find_floridyn_runs(directory::String = pwd())
    if !isdir(directory)
        @warn "Directory does not exist: $directory"
        return String[]
    end
    
    try
        # Get all items in directory
        all_items = readdir(directory)
        
        # Filter for directories starting with "floridyn_run"
        floridyn_dirs = filter(all_items) do item
            full_path = joinpath(directory, item)
            return isdir(full_path) && startswith(item, "floridyn_run")
        end
        
        # Convert to full paths and sort
        full_paths = [joinpath(directory, dir) for dir in floridyn_dirs]
        return sort(full_paths)
        
    catch e
        @error "Error reading directory $directory: $e"
        return String[]
    end
end
