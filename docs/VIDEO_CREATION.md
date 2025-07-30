# Video Creation from Simulation Plots

This document describes how to create videos from saved FLORIDyn simulation plots.

## Overview

FLORIDyn.jl now includes functionality to automatically convert series of PNG plot files into MP4 videos. This is particularly useful for creating animations of flow field evolution over time.

## Prerequisites

- **FFmpeg**: Must be installed and available in your system PATH
  - Ubuntu/Debian: `sudo apt install ffmpeg`
  - MacOS: `brew install ffmpeg`
  - Windows: Download from https://ffmpeg.org/

## Basic Usage

### 1. Generate Plot Frames

First, run your simulation with `vis.save=true` to save plot frames:

```julia
using FLORIDyn

# Enable plot saving
vis = Vis(online=true, save=true, rel_v_min=20.0, up_int=12)

# Run simulation - this will save frames to the 'video/' folder
wf, md, mi = runFLORIDyn(plt, set, wf, wind, sim, con, vis, floridyn, floris)
```

### 2. Create Videos

Create videos from the saved frames:

```julia
# Create video from velocity reduction frames
video_path = createVideo("velocity_reduction"; fps=4)

# Create videos for all measurement types automatically
video_paths = createAllVideos(fps=6, delete_frames=false)

# Custom video creation with options
custom_video = createVideo("wind_speed"; 
                          video_dir="video", 
                          output_dir="animations", 
                          fps=8, 
                          delete_frames=true)
```

## Functions

### `createVideo(prefix; kwargs...)`

Convert PNG files starting with a given prefix into an MP4 video.

**Arguments:**
- `prefix::String`: Prefix of PNG files to include (e.g., "velocity_reduction")
- `video_dir::String`: Input directory (default: "video")
- `output_dir::String`: Output directory (default: "video") 
- `fps::Int`: Frames per second (default: 2)
- `delete_frames::Bool`: Delete PNG files after creating video (default: false)

**Returns:** Path to created video file, or empty string if failed

### `createAllVideos(; kwargs...)`

Automatically create videos for all common measurement types found in the directory.

**Arguments:** Same as `createVideo()` except `prefix`

**Returns:** Vector of paths to created video files

**Supported prefixes:**
- "velocity_reduction" - Velocity reduction animations
- "added_turbulence" - Added turbulence animations  
- "wind_speed" - Wind speed animations

## File Naming Convention

The video creation functions expect PNG files to follow this naming pattern:
- `velocity_reduction_t0000s.png`
- `velocity_reduction_t0012s.png`
- `wind_speed_t0024s.png`
- etc.

Output videos are named:
- `velocity_reduction_animation.mp4`
- `added_turbulence_animation.mp4`
- `wind_speed_animation.mp4`

## Examples

### Example 1: Basic Animation Creation

```julia
using FLORIDyn, ControlPlots

# Run simulation with plot saving
vis = Vis(online=true, save=true)
# ... run simulation ...

# Create animation
video_path = createVideo("velocity_reduction"; fps=4)
println("Created: $video_path")
```

### Example 2: Batch Video Creation

```julia
# Create all available videos
video_paths = createAllVideos(fps=6, delete_frames=false)

for video in video_paths
    println("Created: $video")
end
```

### Example 3: High-Quality Video

```julia
# Create high frame rate video in custom directory
createVideo("wind_speed"; 
           video_dir="simulation_frames",
           output_dir="high_quality_videos", 
           fps=12,
           delete_frames=true)
```

## Tips

1. **Frame Rate**: Start with low FPS (2-4) for overview videos, use higher FPS (8-12) for detailed analysis
2. **Storage**: Use `delete_frames=true` to save disk space after video creation
3. **Quality**: Videos are created with good quality settings (CRF 23) suitable for presentations
4. **Troubleshooting**: If video creation fails, check that FFmpeg is installed and PNG files exist

## Integration with main.jl

The example `main.jl` includes video creation as PLT option 6:

```julia
# Set PLT=6 to create videos from saved frames
PLT = 6
include("examples/main.jl")
```

This will automatically detect and create videos for all available measurement types.
