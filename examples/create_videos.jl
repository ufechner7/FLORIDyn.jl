# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example script demonstrating video creation from simulation plots
using FLORIDyn

# Example 1: Create a video from velocity reduction frames
println("Example 1: Creating video from velocity reduction frames")
video_path = createVideo("velocity_reduction"; fps=6)
if !isempty(video_path)
    println("✓ Created video: $video_path")
else
    println("No velocity reduction frames found or video creation failed")
end

# # Example 2: Create videos for all measurement types
# println("\nExample 2: Creating videos for all measurement types")
# video_paths = createAllVideos(fps=6, delete_frames=false)

# # Example 3: Create video with custom settings
# println("\nExample 3: Custom video creation")
# custom_video = createVideo("wind_speed"; 
#                           video_dir="video", 
#                           output_dir="animations", 
#                           fps=8, 
#                           delete_frames=false)

if !isempty(video_path)
    println("✓ Created custom video: $video_path")
end

println("\nVideo creation examples completed!")
println("Note: FFmpeg must be installed for video creation to work.")
