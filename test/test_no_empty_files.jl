# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Test

@testset "No Empty Files" begin
    @testset "Check for empty files in project" begin
        # Start from project root directory
        project_root = dirname(@__DIR__)  # Parent of test directory
        
        # File extensions to check
        extensions_to_check = [".jl", ".md", ".yml", ".yaml", ".toml"]
        
        # Files that are allowed to be empty (if any)
        allowed_empty_files = [
            # Add any files that are legitimately empty placeholders
            # "src/some_placeholder.jl",
        ]
        
        # Directories to skip (can include hidden dirs that should be ignored)
        skip_directories = [
            ".git",           # Git metadata
            ".github/workflows", # CI files might have empty templates
            "node_modules",   # If any exist
            "coverage",       # Coverage reports might have empty files
            ".vscode",        # VSCode settings
            ".julia",         # Julia environment
        ]
        
        empty_files = String[]
        
        # Function to recursively find all files in a directory (including hidden)
        function find_files_recursive(dir::String, extensions::Vector{String}, skip_dirs::Vector{String})
            files = String[]
            if !isdir(dir)
                return files
            end
            
            try
                for item in readdir(dir; join=false, sort=true)
                    item_path = joinpath(dir, item)
                    
                    # Get relative path from project root for skip checking
                    rel_path = relpath(item_path, project_root)
                    
                    # Skip if this directory is in the skip list
                    if any(startswith(rel_path, skip_dir) for skip_dir in skip_dirs)
                        continue
                    end
                    
                    if isdir(item_path)
                        # Recursively search subdirectories (including hidden ones)
                        append!(files, find_files_recursive(item_path, extensions, skip_dirs))
                    elseif isfile(item_path)
                        # Check if file has one of the target extensions
                        if any(endswith(item_path, ext) for ext in extensions)
                            push!(files, item_path)
                        end
                    end
                end
            catch e
                # Skip directories we can't read (permissions, etc.)
                @warn "Could not read directory $dir: $e"
            end
            return files
        end
        
        # Find all files recursively from project root
        all_files = find_files_recursive(project_root, extensions_to_check, skip_directories)
        
        # Check each file for emptiness
        for file in all_files
            # Convert to relative path for allowed_empty_files comparison
            rel_file = relpath(file, project_root)
            
            # Skip files that are allowed to be empty
            if rel_file in allowed_empty_files
                continue
            end
            
            # Check if file is empty
            try
                if filesize(file) == 0
                    push!(empty_files, rel_file)  # Store relative path for cleaner output
                end
            catch e
                @warn "Could not check file size for $file: $e"
            end
        end
        
        # Test that no empty files were found
        if !isempty(empty_files)
            error_msg = "Found empty files in the project:\n$(join(empty_files, "\n"))"
            println(error_msg)
            @test false
        else
            @test true  # Pass the test if no empty files found
        end
        
        # Print summary
        if isempty(empty_files)
            println("✓ No empty files found in project")
        else
            println("✗ Found $(length(empty_files)) empty file(s)")
        end
    end
end
