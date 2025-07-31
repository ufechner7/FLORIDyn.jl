# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
Video Player Script for FLORIDyn.jl

This script provides functionality to play MP4 videos created by FLORIDyn.jl 
using external video players commonly available on Linux systems.

Usage:
    julia examples/play_video.jl                    # Interactive menu
    julia examples/play_video.jl video_name.mp4     # Play specific video
    julia examples/play_video.jl --list            # List available videos
    julia examples/play_video.jl --all             # Play all videos sequentially
"""

using Printf

# Default video directory where FLORIDyn saves videos
const DEFAULT_VIDEO_DIR = "video"

"""
    detect_video_player()

Detect available video players on the system.
Returns the command for the best available player.
"""
function detect_video_player()
    # List of video players to try, in order of preference
    players = [
        ("mpv", "mpv"),           # Lightweight, high quality
        ("vlc", "vlc"),           # Popular cross-platform player
        ("mplayer", "mplayer"),   # Classic player
        ("totem", "totem"),       # GNOME default
        ("kde-open", "kde-open"), # KDE default
        ("xdg-open", "xdg-open")  # System default
    ]
    
    for (cmd, display_name) in players
        try
            # Check if command exists and is executable
            result = run(pipeline(`which $cmd`, devnull), wait=false)
            wait(result)
            if result.exitcode == 0
                println("🎬 Using video player: $display_name")
                return cmd
            end
        catch
            # Command not found, continue to next
        end
    end
    
    error("❌ No suitable video player found. Please install one of: mpv, vlc, mplayer, totem")
end

"""
    list_videos(video_dir::String=DEFAULT_VIDEO_DIR)

List all MP4 videos in the specified directory.
Returns a vector of video filenames.
"""
function list_videos(video_dir::String=DEFAULT_VIDEO_DIR)
    if !isdir(video_dir)
        @warn "Video directory '$video_dir' does not exist"
        return String[]
    end
    
    videos = filter(f -> endswith(lowercase(f), ".mp4"), readdir(video_dir))
    sort!(videos)  # Sort alphabetically
    return videos
end

"""
    play_video(video_path::String, player_cmd::String)

Play a single video using the specified player command.
"""
function play_video(video_path::String, player_cmd::String)
    if !isfile(video_path)
        @error "Video file not found: $video_path"
        return false
    end
    
    println("▶️  Playing: $(basename(video_path))")
    
    try
        # Run the video player
        if player_cmd == "xdg-open" || player_cmd == "kde-open"
            # For system default players, just open the file
            run(`$player_cmd $video_path`)
        else
            # For dedicated video players, add useful options
            if player_cmd == "mpv"
                # MPV with some nice options
                run(`mpv --loop=inf --osd-fractions --osd-level=2 $video_path`)
            elseif player_cmd == "vlc"
                # VLC without repeat option
                run(pipeline(`vlc --intf dummy --quiet --play-and-exit $video_path`, stdout=devnull, stderr=devnull))
            else
                # Generic player
                run(`$player_cmd $video_path`)
            end
        end
        return true
    catch e
        @error "Failed to play video: $e"
        return false
    end
end

"""
    interactive_menu(video_dir::String=DEFAULT_VIDEO_DIR)

Display an interactive menu to select and play videos.
"""
function interactive_menu(video_dir::String=DEFAULT_VIDEO_DIR)
    videos = list_videos(video_dir)
    
    if isempty(videos)
        println("❌ No MP4 videos found in '$video_dir' directory")
        println("💡 Run your FLORIDyn simulation with video creation enabled first:")
        println("   - Set vis.save=true in your simulation")
        println("   - Use PLT=7 in main.jl to create videos")
        println("   - Or run createVideo() or createAllVideos() functions")
        return
    end
    
    player_cmd = detect_video_player()
    
    while true
        println("\n" * "="^60)
        println("🎬 FLORIDyn Video Player")
        println("="^60)
        println("Available videos in '$video_dir/':")
        println()
        
        for (i, video) in enumerate(videos)
            file_path = joinpath(video_dir, video)
            # Get file size for display
            try
                size_mb = round(stat(file_path).size / 1024 / 1024, digits=1)
                println("  $i. $video ($(size_mb) MB)")
            catch
                println("  $i. $video")
            end
        end
        
        println()
        println("Options:")
        println("  1-$(length(videos)): Play video")
        println("  a: Play all videos sequentially") 
        println("  r: Refresh video list")
        println("  q: Quit")
        println()
        
        print("Enter your choice: ")
        choice = strip(readline())
        
        if choice == "q" || choice == "quit"
            println("👋 Goodbye!")
            break
        elseif choice == "r" || choice == "refresh"
            videos = list_videos(video_dir)
            continue
        elseif choice == "a" || choice == "all"
            println("\n🎬 Playing all videos sequentially...")
            for video in videos
                video_path = joinpath(video_dir, video)
                if !play_video(video_path, player_cmd)
                    break
                end
                
                # Brief pause between videos
                println("\nPress Enter to continue to next video (or Ctrl+C to stop)...")
                try
                    readline()
                catch InterruptException
                    println("\n⏹️  Stopped by user")
                    break
                end
            end
        else
            # Try to parse as video number
            try
                video_num = parse(Int, choice)
                if 1 <= video_num <= length(videos)
                    video_path = joinpath(video_dir, videos[video_num])
                    play_video(video_path, player_cmd)
                else
                    println("❌ Invalid video number. Please choose 1-$(length(videos))")
                end
            catch ArgumentError
                println("❌ Invalid choice. Please enter a number, 'a', 'r', or 'q'")
            end
        end
    end
end

"""
    main()

Main function that handles command line arguments and starts the appropriate functionality.
"""
function main()
    args = ARGS
    
    if isempty(args)
        # No arguments - start interactive menu
        interactive_menu()
    elseif length(args) == 1
        arg = args[1]
        
        if arg == "--help" || arg == "-h"
            println("FLORIDyn Video Player")
            println()
            println("Usage:")
            println("  julia examples/play_video.jl                    # Interactive menu") 
            println("  julia examples/play_video.jl video_name.mp4     # Play specific video")
            println("  julia examples/play_video.jl --list            # List available videos")
            println("  julia examples/play_video.jl --all             # Play all videos")
            println("  julia examples/play_video.jl --help            # Show this help")
            
        elseif arg == "--list" || arg == "-l"
            videos = list_videos()
            if isempty(videos)
                println("❌ No MP4 videos found in '$DEFAULT_VIDEO_DIR' directory")
            else
                println("📁 Available videos in '$DEFAULT_VIDEO_DIR/':")
                for video in videos
                    file_path = joinpath(DEFAULT_VIDEO_DIR, video)
                    try
                        size_mb = round(stat(file_path).size / 1024 / 1024, digits=1)
                        println("  • $video ($(size_mb) MB)")
                    catch
                        println("  • $video")
                    end
                end
            end
            
        elseif arg == "--all" || arg == "-a"
            videos = list_videos()
            if isempty(videos)
                println("❌ No MP4 videos found in '$DEFAULT_VIDEO_DIR' directory")
            else
                player_cmd = detect_video_player()
                println("🎬 Playing all videos sequentially...")
                for video in videos
                    video_path = joinpath(DEFAULT_VIDEO_DIR, video)
                    if !play_video(video_path, player_cmd)
                        break
                    end
                end
            end
            
        else
            # Assume it's a video filename
            video_path = if isfile(arg)
                arg  # Full path provided
            else
                joinpath(DEFAULT_VIDEO_DIR, arg)  # Filename only, look in default directory
            end
            
            if isfile(video_path)
                player_cmd = detect_video_player()
                play_video(video_path, player_cmd)
            else
                @error "Video file not found: $video_path"
                println("💡 Available videos:")
                videos = list_videos()
                for video in videos
                    println("  • $video")
                end
            end
        end
    else
        println("❌ Too many arguments. Use --help for usage information.")
    end
end

main()
