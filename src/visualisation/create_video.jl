# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    createVideo(prefix::String; video_dir="video", output_dir="video", fps=2, delete_frames=false)

Convert PNG files in a directory starting with a given prefix into an MP4 video.

# Arguments
- `prefix::String`: The prefix string that PNG files must start with (e.g., "velocity_reduction", "wind_speed")
- `video_dir::String`: Directory containing the PNG files (default: "video")
- `output_dir::String`: Directory where the output video will be saved (default: "video")
- `fps::Int`: Frames per second for the output video (default: 2)
- `delete_frames::Bool`: Whether to delete the PNG files after creating the video (default: false)

# Returns
- `String`: Path to the created video file, or empty string if creation failed

# Description
This function searches for PNG files in the specified directory that start with the given prefix,
sorts them naturally (handling numeric sequences correctly), and combines them into an MP4 video
using FFmpeg. The function requires FFmpeg to be installed on the system.

# Examples
```julia
# Create video from velocity reduction frames
video_path = createVideo("velocity_reduction"; fps=4)

# Create video from wind speed frames and delete source frames
video_path = createVideo("wind_speed"; fps=6, delete_frames=true)

# Create video from custom directory
video_path = createVideo("added_turbulence"; video_dir="custom_plots", output_dir="videos")
```

# Requirements
- FFmpeg must be installed and available in the system PATH
- PNG files should follow a consistent naming pattern with the prefix
- Recommended naming: "prefix_t0000s.png", "prefix_t0012s.png", etc.

# Notes
- Files are sorted naturally to handle numeric sequences correctly (e.g., t0001s, t0010s, t0100s)
- The output video filename will be "prefix_animation.mp4"
- If no matching files are found, the function returns an empty string
- FFmpeg parameters are optimized for good quality and reasonable file size
"""
function createVideo(prefix::String; video_dir="video", output_dir="video", fps=2, delete_frames=false)
    # Check if video directory exists
    if !isdir(video_dir)
        @warn "Video directory '$video_dir' does not exist"
        return ""
    end
    
    # Find all PNG files starting with the prefix
    all_files = readdir(video_dir)
    matching_files = filter(f -> startswith(f, prefix) && endswith(f, ".png"), all_files)
    
    if isempty(matching_files)
        @warn "No PNG files found with prefix '$prefix' in directory '$video_dir'"
        return ""
    end
    
    # Sort files naturally (handles numeric sequences correctly)
    sorted_files = sort(matching_files, by=natural_sort_key)
    
    println("Found $(length(sorted_files)) files with prefix '$prefix'")
    println("First file: $(sorted_files[1])")
    println("Last file: $(sorted_files[end])")
    
    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # Generate output filename
    output_filename = "$(prefix)_animation.mp4"
    output_path = joinpath(output_dir, output_filename)
    
    # Create temporary file list for FFmpeg
    temp_list_file = "temp_file_list.txt"
    
    try
        # Write file list for FFmpeg concat demuxer
        open(temp_list_file, "w") do f
            for file in sorted_files
                full_path = joinpath(video_dir, file)
                # FFmpeg requires forward slashes even on Windows
                normalized_path = replace(full_path, "\\" => "/")
                println(f, "file '$normalized_path'")
                println(f, "duration $(1.0/fps)")
            end
            # Add the last frame one more time to ensure it's displayed
            if !isempty(sorted_files)
                last_file = joinpath(video_dir, sorted_files[end])
                normalized_path = replace(last_file, "\\" => "/")
                println(f, "file '$normalized_path'")
            end
        end
        
        # Build FFmpeg command
        ffmpeg_cmd = `ffmpeg -y -f concat -safe 0 -i $temp_list_file -vf "fps=$fps,scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -preset medium -crf 23 -pix_fmt yuv420p $output_path`
        
        println("Creating video with FFmpeg...")
        println("Command: $ffmpeg_cmd")
        
        # Run FFmpeg
        try
            run(ffmpeg_cmd)
            println("✓ Video created successfully: $output_path")
            
            # Delete PNG files if requested
            if delete_frames
                println("Deleting source PNG files...")
                for file in sorted_files
                    file_path = joinpath(video_dir, file)
                    try
                        rm(file_path)
                    catch e
                        @warn "Failed to delete file $file_path: $e"
                    end
                end
                println("✓ Source files deleted")
            end
            
            return output_path
            
        catch e
            if isa(e, ProcessFailedException)
                @error "FFmpeg failed to create video. Make sure FFmpeg is installed and available in PATH."
                @error "Error details: $e"
            else
                @error "Unexpected error running FFmpeg: $e"
            end
            return ""
        end
        
    catch e
        @error "Error creating video: $e"
        return ""
        
    finally
        # Clean up temporary file
        if isfile(temp_list_file)
            try
                rm(temp_list_file)
            catch
                # Ignore cleanup errors
            end
        end
    end
end

"""
    natural_sort_key(filename::String)

Generate a sort key for natural sorting of filenames containing numbers.

# Arguments
- `filename::String`: The filename to generate a sort key for

# Returns
- `Vector`: Sort key that handles numeric sequences naturally

# Description
This function creates a sort key that handles numeric sequences in filenames correctly.
For example, it will sort ["file1.png", "file10.png", "file2.png"] as 
["file1.png", "file2.png", "file10.png"] rather than alphabetically.

# Examples
```julia
files = ["velocity_reduction_t0001s.png", "velocity_reduction_t0010s.png", "velocity_reduction_t0002s.png"]
sorted_files = sort(files, by=natural_sort_key)
# Result: ["velocity_reduction_t0001s.png", "velocity_reduction_t0002s.png", "velocity_reduction_t0010s.png"]
```
"""
function natural_sort_key(filename::String)
    # Split filename into alternating text and number parts
    parts = []
    current_part = ""
    in_number = false
    
    for char in filename
        if isdigit(char)
            if !in_number
                # Starting a number sequence, save previous text part
                if !isempty(current_part)
                    push!(parts, current_part)
                end
                current_part = string(char)
                in_number = true
            else
                # Continue number sequence
                current_part *= char
            end
        else
            if in_number
                # Ending a number sequence, save as integer
                if !isempty(current_part)
                    push!(parts, parse(Int, current_part))
                end
                current_part = string(char)
                in_number = false
            else
                # Continue text sequence
                current_part *= char
            end
        end
    end
    
    # Add the final part
    if !isempty(current_part)
        if in_number
            push!(parts, parse(Int, current_part))
        else
            push!(parts, current_part)
        end
    end
    
    return parts
end

"""
    createAllVideos(; video_dir="video", output_dir="video", fps=2, delete_frames=false)

Create videos for all common measurement types found in the video directory.

# Arguments
- `video_dir::String`: Directory containing the PNG files (default: "video")
- `output_dir::String`: Directory where output videos will be saved (default: "video")
- `fps::Int`: Frames per second for output videos (default: 2)
- `delete_frames::Bool`: Whether to delete PNG files after creating videos (default: false)

# Returns
- `Vector{String}`: Paths to created video files

# Description
This convenience function automatically detects common measurement type prefixes in the video
directory and creates videos for each type found. It looks for the following prefixes:
- "velocity_reduction"
- "added_turbulence" 
- "wind_speed"

# Example
```julia
# Create videos for all measurement types found
video_paths = createAllVideos(fps=4, delete_frames=true)
println("Created videos: ", video_paths)
```
"""
function createAllVideos(; video_dir="video", output_dir="video", fps=2, delete_frames=false)
    # Common measurement type prefixes
    prefixes = ["velocity_reduction", "added_turbulence", "wind_speed"]
    created_videos = String[]
    
    println("Searching for videos to create...")
    
    for prefix in prefixes
        video_path = createVideo(prefix; video_dir=video_dir, output_dir=output_dir, fps=fps, delete_frames=delete_frames)
        if !isempty(video_path)
            push!(created_videos, video_path)
        end
    end
    
    if isempty(created_videos)
        println("No videos were created - no matching PNG files found")
    else
        println("\n✓ Successfully created $(length(created_videos)) video(s):")
        for video in created_videos
            println("  - $video")
        end
    end
    
    return created_videos
end
