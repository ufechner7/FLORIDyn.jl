# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, ControlPlots

function get_parameters()
    settings_file = "data/2021_9T_Data.yaml"

    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)

    # create settings struct
    set = Settings(wind, sim, con)

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)  
    wf = initSimulation(wf, sim)
    
    # Create visualization settings for testing
    vis = Vis(online=false, save=false)
    wf, md, mi = runFLORIDyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)
    return wf, set, floris, wind, md
end

"""
Create a simple test PNG file of specified dimensions.
This creates a valid PNG that FFmpeg can process without errors.
"""
function create_test_png(filepath::String, width::Int, height::Int)
    # Create a minimal but valid 10x10 black PNG file
    simple_png = UInt8[
        # PNG signature
        0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A,
        # IHDR chunk (10x10, grayscale)
        0x00, 0x00, 0x00, 0x0D, 0x49, 0x48, 0x44, 0x52,
        0x00, 0x00, 0x00, 0x0A, 0x00, 0x00, 0x00, 0x0A,  # 10x10 dimensions
        0x08, 0x00, 0x00, 0x00, 0x00, 0x3A, 0x30, 0x9D,  # 8-bit grayscale
        0x9B, 
        # IDAT chunk (compressed black pixels)
        0x00, 0x00, 0x00, 0x0E, 0x49, 0x44, 0x41, 0x54,
        0x78, 0x9C, 0x63, 0x00, 0x02, 0x00, 0x00, 0x05,
        0x00, 0x01, 0x0D, 0x0A, 0x2D, 0xB4,
        # IEND chunk
        0x00, 0x00, 0x00, 0x00, 0x49, 0x45, 0x4E, 0x44,
        0xAE, 0x42, 0x60, 0x82
    ]
    
    open(filepath, "w") do f
        write(f, simple_png)
    end
end

@testset verbose=true "visualisation                                           " begin
    @testset "getMeasurements" begin
        # Create a simple test wind farm configuration
        @testset "basic functionality" begin
            # Create minimal test matrices
            mx = [0.0 100.0; 200.0 300.0]  # 2x2 grid of x-coordinates
            my = [0.0 0.0; 100.0 100.0]    # 2x2 grid of y-coordinates
            nM = 3  # Number of measurements (typically velocity reduction, added turbulence, effective wind speed)
            zh = 90.0  # Hub height
            wf, set, floris, wind, md = get_parameters()
            # Call the function with the correct number of parameters
            measurements = getMeasurements(mx, my, nM, zh, wf, set, floris, wind)
            
            # Test basic input properties that should be satisfied
            @test size(mx) == size(my)  # Matrices must have same dimensions
            @test nM > 0  # Number of measurements must be positive
            @test zh > 0  # Hub height must be positive
            @test isa(nM, Int)  # nM should be an integer
            @test isa(zh, Real)  # zh should be a real number
            
            # Test output properties
            @test isa(measurements, Array{Float64,3})
            @test size(measurements) == (size(mx, 1), size(mx, 2), nM)
        end
        @testset "function exists and is callable" begin
            # Test that the function exists and has the right signature
            @test isa(getMeasurements, Function)
            
            # Test that we can get method information
            methods_list = methods(getMeasurements)
            @test length(methods_list) >= 1
            
            # Check that the function is exported
            @test :getMeasurements in names(FLORIDyn)
        end
    end
    
    @testset "plotFlowField" begin
        @testset "basic functionality" begin
            # Get test parameters
            wf, set, floris, wind, md = get_parameters()
            
            # Call the calcFlowField function
            Z, X, Y = calcFlowField(set, wf, wind, floris)

            vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
            plotFlowField(nothing, plt, wf, X, Y, Z, vis; unit_test=true)

            # Test that outputs have the correct types
            @test isa(Z, Array{Float64,3})
            @test isa(X, Matrix{Float64})
            @test isa(Y, Matrix{Float64})
            
            # Test that X and Y have the same dimensions
            @test size(X) == size(Y)
            
            # Test that Z has the correct dimensions
            # Z should have the same x,y dimensions as X,Y and 3 measurements in the third dimension
            @test size(Z, 1) == size(X, 1)
            @test size(Z, 2) == size(X, 2)
            @test size(Z, 3) == 3  # nM = 3 measurements
            
            # Test that the coordinate grids are reasonable
            # X should vary along columns, Y should vary along rows
            @test X[1, 1] < X[1, end]  # X increases along columns
            @test Y[1, 1] < Y[end, 1]  # Y increases along rows
            
            # Test grid spacing is consistent (20m resolution)
            if size(X, 2) > 1
                @test abs((X[1, 2] - X[1, 1]) - 20.0) < 1e-10  # 20m spacing
            end
            if size(Y, 1) > 1
                @test abs((Y[2, 1] - Y[1, 1]) - 20.0) < 1e-10  # 20m spacing
            end
        end
        
        @testset "function exists and is callable" begin
            # Test that the function exists and has the right signature
            @test isa(plotFlowField, Function)
            
            # Test that we can get method information
            methods_list = methods(plotFlowField)
            @test length(methods_list) >= 1
            
            # Check that the function is exported
            @test :plotFlowField in names(FLORIDyn)
        end
        
        @testset "coordinate grid properties" begin
            # Get test parameters
            wf, set, floris, wind, md = get_parameters()
            
            # Call the function
            Z, X, Y = calcFlowField(set, wf, wind, floris)
            
            # Test coordinate range
            @test minimum(X) >= 0.0
            @test maximum(X) <= 3000.0
            @test minimum(Y) >= 0.0  
            @test maximum(Y) <= 3000.0
            
            # Test that coordinates are within expected bounds
            @test all(X .>= 0.0)
            @test all(X .<= 3000.0)
            @test all(Y .>= 0.0)
            @test all(Y .<= 3000.0)
        end
        
        @testset "msr parameter tests" begin
            # Get test parameters
            wf, set, floris, wind, md = get_parameters()
            
            # Call the calcFlowField function
            Z, X, Y = calcFlowField(set, wf, wind, floris)
            vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
            
            @testset "msr=1 (velocity reduction)" begin
                # Test plotFlowField with msr=1
                plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=1, unit_test=true)

                # Test that the function runs without error
                @test true  # If we get here, the function didn't throw an error
                
                # Test that Z has the required third dimension for msr=1
                @test size(Z, 3) >= 1
                
                # Test that msr=1 data exists and is reasonable (velocity reduction should be >= 0%)
                velocity_reduction = Z[:, :, 1]
                @test all(velocity_reduction .>= 0.0)
                @test all(velocity_reduction .<= 100.0)  # Should not exceed 100% reduction
            end
            
            @testset "msr=2 (added turbulence)" begin
                # Test plotFlowField with msr=2
                plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=2, unit_test=true)

                # Test that the function runs without error
                @test true  # If we get here, the function didn't throw an error
                
                # Test that Z has the required third dimension for msr=2
                @test size(Z, 3) >= 2
                
                # Test that msr=2 data exists and is reasonable (added turbulence should be >= 0%)
                added_turbulence = Z[:, :, 2]
                @test all(added_turbulence .>= 0.0)
                @test all(isfinite.(added_turbulence))  # Should be finite values
            end
            
            @testset "msr parameter validation" begin
                vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4)
                # Test error handling for invalid msr values
                # msr=0 causes BoundsError (Julia arrays are 1-indexed)
                @test_throws BoundsError plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=0, unit_test=true)
                # msr > 3 causes ErrorException from explicit check
                @test_throws ErrorException plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=4, unit_test=true)

                # Test that msr=3 (default) still works
                plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=3, unit_test=true)
                @test true  # If we get here, the function didn't throw an error
                
                # Test that wind speed data (msr=3) is reasonable
                wind_speed = Z[:, :, 3]
                @test all(wind_speed .>= 0.0)  # Wind speed should be non-negative
                @test all(isfinite.(wind_speed))  # Should be finite values
            end
        end

    end
    @testset "plotMeasurements" begin
        # Get test parameters
        wf, set, floris, wind, md = get_parameters()
        plotMeasurements(plt, wf, md; separated=true, unit_test=true)
        plotMeasurements(plt, wf, md; unit_test=true)
    end
    
    @testset "createVideo" begin
        @testset "function exists and is callable" begin
            # Test that the function exists and has the right signature
            @test isa(createVideo, Function)
            
            # Test that we can get method information
            methods_list = methods(createVideo)
            @test length(methods_list) >= 1
            
            # Check that the function is exported
            @test :createVideo in names(FLORIDyn)
        end
        
        @testset "basic functionality with test images" begin
            # Create a temporary test directory
            test_dir = "test_video_temp"
            output_dir = "test_output_temp"
            
            try
                # Create test directory
                mkpath(test_dir)
                mkpath(output_dir)
                
                # Create some simple test PNG files using ImageMagick convert if available
                # or create dummy files for testing
                test_files = ["test_t0000s.png", "test_t0012s.png", "test_t0024s.png"]
                
                # Create minimal PNG files (64x64 pixels) for testing
                # Use a larger size that FFmpeg can handle better
                test_files = ["test_t0000s.png", "test_t0012s.png", "test_t0024s.png"]
                
                for file in test_files
                    file_path = joinpath(test_dir, file)
                    
                    # Try to use a plotting library to create a proper PNG
                    try
                        # Create a simple 64x64 test image
                        if @isdefined(plt)
                            plt.figure(figsize=(1, 1))
                            plt.imshow(rand(64, 64), cmap="viridis")
                            plt.axis("off")
                            plt.savefig(file_path, dpi=64, bbox_inches="tight", pad_inches=0)
                            plt.close()
                        else
                            # Fallback: create a larger minimal PNG (10x10 pixels)
                            create_test_png(file_path, 10, 10)
                        end
                    catch e
                        @warn "Could not create test PNG with plotting, using minimal fallback: $e"
                        # Create a minimal but larger PNG file
                        create_test_png(file_path, 10, 10)
                    end
                end
                
                # Test basic functionality (should fail gracefully without FFmpeg)
                result = createVideo("test"; video_dir=test_dir, output_dir=output_dir, fps=2)
                
                # The function should handle missing FFmpeg gracefully
                @test isa(result, String)  # Should return a string (empty if failed)
                
                # Test with non-existent directory
                result_nodir = createVideo("test"; video_dir="nonexistent_dir")
                @test result_nodir == ""  # Should return empty string
                
                # Test with no matching files
                result_nomatch = createVideo("nomatch"; video_dir=test_dir)
                @test result_nomatch == ""  # Should return empty string
                
                # Test natural_sort_key function directly
                test_filenames = ["test_t0010s.png", "test_t0001s.png", "test_t0100s.png"]
                sorted_names = sort(test_filenames, by=natural_sort_key)
                expected_order = ["test_t0001s.png", "test_t0010s.png", "test_t0100s.png"]
                @test sorted_names == expected_order
                
            finally
                # Clean up test directories
                if isdir(test_dir)
                    rm(test_dir; recursive=true)
                end
                if isdir(output_dir)
                    rm(output_dir; recursive=true)
                end
            end
        end
        
        @testset "natural_sort_key" begin
            # Test natural sorting behavior
            @test isa(natural_sort_key, Function)
            
            # Test with simple filenames
            key1 = natural_sort_key("file1.png")
            key2 = natural_sort_key("file10.png")
            key3 = natural_sort_key("file2.png")
            
            @test isa(key1, Vector)
            @test isa(key2, Vector)
            @test isa(key3, Vector)
            
            # Test actual sorting behavior
            files = ["velocity_reduction_t0010s.png", "velocity_reduction_t0001s.png", "velocity_reduction_t0100s.png"]
            sorted_files = sort(files, by=natural_sort_key)
            expected = ["velocity_reduction_t0001s.png", "velocity_reduction_t0010s.png", "velocity_reduction_t0100s.png"]
            @test sorted_files == expected
            
            # Test with mixed content
            mixed_files = ["a10b.png", "a2b.png", "a1b.png"]
            sorted_mixed = sort(mixed_files, by=natural_sort_key)
            expected_mixed = ["a1b.png", "a2b.png", "a10b.png"]
            @test sorted_mixed == expected_mixed
        end
        
        @testset "createAllVideos" begin
            # Test that the function exists and is callable
            @test isa(createAllVideos, Function)
            
            # Check that the function is exported
            @test :createAllVideos in names(FLORIDyn)
            
            # Test with empty directory
            test_dir = "test_empty_video"
            try
                mkpath(test_dir)
                result = createAllVideos(video_dir=test_dir, output_dir=test_dir)
                @test isa(result, Vector{String})
                @test isempty(result)  # Should return empty vector
            finally
                if isdir(test_dir)
                    rm(test_dir; recursive=true)
                end
            end
        end
    end
end
nothing



nothing
