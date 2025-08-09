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

Delete the last (newest) n folders starting with "floridyn_run" from the output directory.

This function helps manage disk space by removing recent simulation run directories. It sorts
directories by modification time and removes the newest ones first, preserving older
simulation results that may be more established or important. The function operates on the
parent directory of the visualization output path.

# Arguments
- `vis::Vis`: Visualization settings object containing output path information
- `n::Int`: Number of folders to delete (must be positive, default: 1)
- `dry_run::Bool`: If `true`, only shows what would be deleted without actually deleting (default: `false`)

# Returns
- `Vector{String}`: List of directories that were deleted (or would be deleted in dry_run mode)

# Examples
```julia
# Create visualization settings
vis = Vis("data/vis_default.yaml")

# Delete the newest floridyn_run folder from vis.output_path directory
deleted = delete_results(vis)

# Delete the 3 newest floridyn_run folders
deleted = delete_results(vis, 3)

# Preview what would be deleted without actually deleting
delete_results(vis, 5, dry_run=true)
```

# Behavior
- Automatically determines the directory to search from `dirname(vis.output_path)`
- Sets `vis.unique_folder = ""` as part of the cleanup process
- Only deletes directories that start with "floridyn_run"
- Sorts by modification time (newest deleted first)

# Error Handling
- Returns empty list if no matching directories found
- Skips directories that cannot be deleted due to permissions
- Validates positive number of folders to delete

# Notes
- Preserves the oldest simulation results
- Useful for removing failed or incomplete recent runs, or test runs
- Works with directories created by [`unique_name()`](@ref)
- Integrates with the [`Vis`](@ref) settings system

# See Also
- [`Vis`](@ref): Visualization settings struct
- [`unique_name`](@ref): Function that creates floridyn_run directories
- [`find_floridyn_runs`](@ref): Function to list existing floridyn_run directories
"""
function delete_results(vis::Vis, n::Int=1; dry_run::Bool = false)
    # Validate input
    vis.unique_folder = ""
    directory = vis.output_path
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
    
    # Actually delete the directories
    deleted_dirs = String[]
    
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
