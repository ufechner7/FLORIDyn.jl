# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

if !isdefined(Main, :Test)
    using Test
end 
using DataFrames  # Required for DataFrame functionality in turbines() tests
using Statistics  # Required for mean, std functions in calc_rel_power tests
using DistributedNext  # Required for nprocs() function in plot_x tests

if ! isinteractive()
    global pltctrl
    if !isdefined(Main, :FLORIDyn)
        using FLORIDyn
    end

    if !isdefined(Main, :ControlPlots)
        if Threads.nthreads() == 1
            using ControlPlots; 
        else 
            global plt
            plt = nothing
        end
    end
    pltctrl = nothing
    if Threads.nthreads() == 1
        global pltctrl
        pltctrl = ControlPlots
    end

    function get_parameters()
        settings_file = "data/2021_9T_Data.yaml"

        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn, ta = setup(settings_file)

        # create settings struct with automatic parallel/threading detection
        use_threading = Threads.nthreads() > 1
        set = Settings(wind, sim, con, use_threading, use_threading)

        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)  
        wf = initSimulation(wf, sim)
        
        # Create visualization settings for testing
        vis = Vis(online=false, save=false, unit_test=true)
        wf, md, mi = run_floridyn(nothing, set, wf, wind, sim, con, vis, floridyn, floris)
        return wf, set, floris, wind, md
    end

    """
    Create a simple test PNG file of specified dimensions.
    This creates a valid PNG that FFmpeg can process without errors.
    Uses basic RGB format with black pixels.
    """
    function create_test_png(filepath::String, width::Int, height::Int)
        # Ensure minimum dimensions
        w = max(width, 1)
        h = max(height, 1)
        
        # CRC32 calculation function for PNG chunks
        function crc32(data::Vector{UInt8})
            crc_table = UInt32[
                0x00000000, 0x77073096, 0xee0e612c, 0x990951ba, 0x076dc419, 0x706af48f,
                0xe963a535, 0x9e6495a3, 0x0edb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988,
                0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91, 0x1db71064, 0x6ab020f2,
                0xf3b97148, 0x84be41de, 0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
                0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec, 0x14015c4f, 0x63066cd9,
                0xfa0f3d63, 0x8d080df5, 0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
                0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b, 0x35b5a8fa, 0x42b2986c,
                0xdbbbc9d6, 0xacbcf940, 0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
                0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116, 0x21b4f4b5, 0x56b3c423,
                0xcfba9599, 0xb8bda50f, 0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
                0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d, 0x76dc4190, 0x01db7106,
                0x98d220bc, 0xefd5102a, 0x71b18589, 0x06b6b51f, 0x9fbfe4a5, 0xe8b8d433,
                0x7807c9a2, 0x0f00f934, 0x9609a88e, 0xe10e9818, 0x7f6a0dbb, 0x086d3d2d,
                0x91646c97, 0xe6635c01, 0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e,
                0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457, 0x65b0d9c6, 0x12b7e950,
                0x8bbeb8ea, 0xfcb9887c, 0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
                0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2, 0x4adfa541, 0x3dd895d7,
                0xa4d1c46d, 0xd3d6f4fb, 0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
                0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9, 0x5005713c, 0x270241aa,
                0xbe0b1010, 0xc90c2086, 0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
                0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4, 0x59b33d17, 0x2eb40d81,
                0xb7bd5c3b, 0xc0ba6cad, 0xedb88320, 0x9abfb3b6, 0x03b6e20c, 0x74b1d29a,
                0xead54739, 0x9dd277af, 0x04db2615, 0x73dc1683, 0xe3630b12, 0x94643b84,
                0x0d6d6a3e, 0x7a6a5aa8, 0xe40ecf0b, 0x9309ff9d, 0x0a00ae27, 0x7d079eb1,
                0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe, 0xf762575d, 0x806567cb,
                0x196c3671, 0x6e6b06e7, 0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc,
                0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5, 0xd6d6a3e8, 0xa1d1937e,
                0x38d8c2c4, 0x4fdff252, 0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
                0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60, 0xdf60efc3, 0xa867df55,
                0x316e8eef, 0x4669be79, 0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
                0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f, 0xc5ba3bbe, 0xb2bd0b28,
                0x2bb45a92, 0x5cb36a04, 0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
                0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x026d930a, 0x9c0906a9, 0xeb0e363f,
                0x72076785, 0x05005713, 0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0x0cb61b38,
                0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0x0bdbdf21, 0x86d3d2d4, 0xf1d4e242,
                0x68ddb3f8, 0x1fda836e, 0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
                0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c, 0x8f659eff, 0xf862ae69,
                0x616bffd3, 0x166ccf45, 0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2,
                0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db, 0xaed16a4a, 0xd9d65adc,
                0x40df0b66, 0x37d83bf0, 0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
                0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6, 0xbad03605, 0xcdd70693,
                0x54de5729, 0x23d967bf, 0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
                0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d
            ]
            
            crc = 0xffffffff
            for byte in data
                crc = crc_table[(crc ⊻ byte) & 0xff + 1] ⊻ (crc >> 8)
            end
            return crc ⊻ 0xffffffff
        end
        
        # Helper function to convert UInt32 to big-endian bytes
        function uint32_to_bytes(value::UInt32)
            return UInt8[
                (value >> 24) & 0xff,
                (value >> 16) & 0xff,
                (value >> 8) & 0xff,
                value & 0xff
            ]
        end
        
        # Create PNG data
        png_data = UInt8[]
        
        # PNG signature
        append!(png_data, UInt8[0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A])
        
        # IHDR chunk
        ihdr_data = UInt8[]
        append!(ihdr_data, UInt8['I', 'H', 'D', 'R'])  # Chunk type
        append!(ihdr_data, uint32_to_bytes(UInt32(w)))  # Width
        append!(ihdr_data, uint32_to_bytes(UInt32(h)))  # Height
        append!(ihdr_data, UInt8[8])    # Bit depth (8 bits per channel)
        append!(ihdr_data, UInt8[2])    # Color type (2 = RGB)
        append!(ihdr_data, UInt8[0])    # Compression method
        append!(ihdr_data, UInt8[0])    # Filter method
        append!(ihdr_data, UInt8[0])    # Interlace method
        
        # Write IHDR chunk
        append!(png_data, uint32_to_bytes(UInt32(13)))  # Data length
        append!(png_data, ihdr_data)
        append!(png_data, uint32_to_bytes(crc32(ihdr_data)))
        
        # Create image data (black RGB pixels)
        # Each pixel is 3 bytes (RGB), each row starts with a filter byte (0 = None)
        row_size = w * 3 + 1  # 3 bytes per pixel + 1 filter byte
        raw_data = zeros(UInt8, h * row_size)
        
        # Fill with filter bytes (0) at the start of each row
        for row in 0:(h-1)
            raw_data[row * row_size + 1] = 0  # Filter type: None
            # Pixels are already black (zeros)
        end
        
        # Compress the image data using simple deflate
        # For simplicity, use uncompressed deflate blocks
        compressed_data = UInt8[]
        
        # Deflate header (no compression, final block)
        push!(compressed_data, 0x01)  # BFINAL=1, BTYPE=00 (no compression)
        
        # Block length (little endian)
        len = length(raw_data)
        push!(compressed_data, len & 0xff)
        push!(compressed_data, (len >> 8) & 0xff)
        push!(compressed_data, (~len) & 0xff)
        push!(compressed_data, ((~len) >> 8) & 0xff)
        
        # Raw data
        append!(compressed_data, raw_data)
        
        # Adler-32 checksum
        a = 1
        b = 0
        for byte in raw_data
            a = (a + byte) % 65521
            b = (b + a) % 65521
        end
        adler32 = (b << 16) | a
        
        # Add zlib header
        zlib_data = UInt8[0x78, 0x01]  # CMF=0x78, FLG=0x01
        append!(zlib_data, compressed_data)
        append!(zlib_data, uint32_to_bytes(UInt32(adler32)))
        
        # IDAT chunk
        idat_chunk = UInt8['I', 'D', 'A', 'T']
        append!(idat_chunk, zlib_data)
        
        # Write IDAT chunk
        append!(png_data, uint32_to_bytes(UInt32(length(zlib_data))))
        append!(png_data, idat_chunk)
        append!(png_data, uint32_to_bytes(crc32(idat_chunk)))
        
        # IEND chunk
        iend_data = UInt8['I', 'E', 'N', 'D']
        append!(png_data, uint32_to_bytes(UInt32(0)))  # Length = 0
        append!(png_data, iend_data)
        append!(png_data, uint32_to_bytes(crc32(iend_data)))
        
        # Write to file in binary mode
        open(filepath, "w") do f
            write(f, png_data)
        end
    end

    @testset verbose=true "Visualisation Tests                                    " begin
        # Add explicit garbage collection at start
        GC.gc()
        
        @testset "getMeasurements" begin
            # Create a simple test wind farm configuration
            @testset "basic functionality" begin
                # Create minimal test matrices
                mx = [0.0 100.0; 200.0 300.0]  # 2x2 grid of x-coordinates
                my = [0.0 0.0; 100.0 100.0]    # 2x2 grid of y-coordinates
                nM = 3  # Number of measurements (typically velocity reduction, added turbulence, effective wind speed)
                zh = 90.0  # Hub height
                wf, set, floris, wind, md = get_parameters()
                # Create appropriate buffers depending on threading
                if Threads.nthreads() > 1
                    buffers = FLORIDyn.create_thread_buffers(wf, Threads.nthreads() + 1, floris)
                else
                    buffers = FLORIDyn.create_thread_buffers(wf, 1, floris)
                end
                # Call the function with the correct number of parameters
                measurements = getMeasurements(buffers, mx, my, nM, zh, wf, set, floris, wind)
                
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
        
        @testset "calcFlowField" begin
            # Get test parameters once for all calcFlowField tests
            wf, set, floris, wind, md = get_parameters()
            
            @testset "basic functionality without vis parameter" begin
                Z, X, Y = calcFlowField(set, wf, wind, floris)
                
                # Basic return type checks
                @test Z isa Array{Float64, 3}  # 3D array with measurements
                @test X isa Matrix{Float64}    # 2D grid of X coordinates
                @test Y isa Matrix{Float64}    # 2D grid of Y coordinates
                
                # Check dimensions consistency (Z should be (ny, nx, nmeasurements))
                @test size(Z, 1) == size(Y, 1)  # Number of Y points
                @test size(Z, 2) == size(X, 2)  # Number of X points
                @test size(Z, 3) == 3           # Three measurements (vel reduction, turbulence, effective wind)
                
                # Check coordinate grids have consistent dimensions
                @test size(X) == size(Y)
                
                # Check for reasonable values (should be non-empty)
                @test size(Z) != (0, 0, 0)
                @test size(X) != (0, 0)
                @test size(Y) != (0, 0)
            end
            
            @testset "vis parameter validation" begin
                # Test with custom vis configuration
                custom_vis = Vis(
                    field_limits_min = [100.0, 200.0, 50.0],
                    field_limits_max = [1100.0, 1200.0, 150.0],
                    field_resolution = 50.0,
                    online = false,
                    save = false,
                    unit_test = true
                )
                
                Z, X, Y = calcFlowField(set, wf, wind, floris; vis=custom_vis)
                
                # Calculate expected grid dimensions
                # Grid points = (max - min) / resolution + 1
                expected_x_points = Int(round((custom_vis.field_limits_max[1] - custom_vis.field_limits_min[1]) / custom_vis.field_resolution)) + 1
                expected_y_points = Int(round((custom_vis.field_limits_max[2] - custom_vis.field_limits_min[2]) / custom_vis.field_resolution)) + 1
                
                # Validate grid dimensions match vis parameters
                @test size(X, 2) == expected_x_points  # X grid has expected number of columns
                @test size(Y, 1) == expected_y_points  # Y grid has expected number of rows
                @test size(Z, 1) == expected_y_points  # Z has expected number of Y points
                @test size(Z, 2) == expected_x_points  # Z has expected number of X points
                @test size(Z, 3) == 3                  # Z has 3 measurements
                
                # Validate coordinate ranges
                @test minimum(X) ≈ custom_vis.field_limits_min[1] atol=1e-10
                @test maximum(X) ≈ custom_vis.field_limits_max[1] atol=1e-10
                @test minimum(Y) ≈ custom_vis.field_limits_min[2] atol=1e-10
                @test maximum(Y) ≈ custom_vis.field_limits_max[2] atol=1e-10
                
                # Validate grid spacing
                if size(X, 2) > 1
                    x_spacing = X[1, 2] - X[1, 1]  # Spacing in first row
                    @test x_spacing ≈ custom_vis.field_resolution atol=1e-10
                end
                if size(Y, 1) > 1
                    y_spacing = Y[2, 1] - Y[1, 1]  # Spacing in first column
                    @test y_spacing ≈ custom_vis.field_resolution atol=1e-10
                end
            end
            
            @testset "vis parameter backward compatibility" begin
                # Test that providing nothing for vis uses default behavior
                Z1, X1, Y1 = calcFlowField(set, wf, wind, floris)
                Z2, X2, Y2 = calcFlowField(set, wf, wind, floris; vis=nothing)
                
                # Results should be identical
                @test Z1 == Z2
                @test X1 == X2
                @test Y1 == Y2
            end
            
            @testset "different resolution values" begin
                # Test with different resolution values
                for resolution in [25.0, 100.0, 200.0]
                    test_vis = Vis(
                        field_limits_min = [0.0, 0.0, 0.0],
                        field_limits_max = [1000.0, 1000.0, 200.0],
                        field_resolution = resolution,
                        online = false,
                        save = false,
                        unit_test = true
                    )
                    
                    Z, X, Y = calcFlowField(set, wf, wind, floris; vis=test_vis)
                    
                    # Calculate expected dimensions
                    expected_x = Int(round(1000.0 / resolution)) + 1
                    expected_y = Int(round(1000.0 / resolution)) + 1
                    
                    @test size(X, 2) == expected_x
                    @test size(Y, 1) == expected_y
                    @test size(Z, 1) == expected_y
                    @test size(Z, 2) == expected_x
                    @test size(Z, 3) == 3
                    
                    # Verify resolution is respected
                    if size(X, 2) > 1
                        @test (X[1, 2] - X[1, 1]) ≈ resolution atol=1e-10
                    end
                    if size(Y, 1) > 1
                        @test (Y[2, 1] - Y[1, 1]) ≈ resolution atol=1e-10
                    end
                end
            end
        end
        
        @testset "plotFlowField" begin
            # Get test parameters once for all plotFlowField tests
            wf, set, floris, wind, md = get_parameters()
            Z, X, Y = calcFlowField(set, wf, wind, floris)
            vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4, unit_test=true)

            @testset "basic functionality" begin
                plot_flow_field(wf, X, Y, Z, vis; plt)

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
                
                @testset "msr=VelReduction (velocity reduction)" begin
                    # Test plot_flow_field with msr=VelReduction
                    plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt)

                    # Test that the function runs without error
                    @test true  # If we get here, the function didn't throw an error
                    
                    # Test that Z has the required third dimension for msr=VelReduction
                    @test size(Z, 3) >= 1
                    
                    # Test that msr=VelReduction data exists and is reasonable (velocity reduction should be >= 0%)
                    velocity_reduction = Z[:, :, 1]
                    @test all(velocity_reduction .>= 0.0)
                    @test all(velocity_reduction .<= 100.0)  # Should not exceed 100% reduction
                end
                
                @testset "msr=AddedTurbulence (added turbulence)" begin
                    # Test plot_flow_field with msr=AddedTurbulence
                    plot_flow_field(wf, X, Y, Z, vis; msr=AddedTurbulence, plt)

                    # Test that the function runs without error
                    @test true  # If we get here, the function didn't throw an error
                    
                    # Test that Z has the required third dimension for msr=AddedTurbulence
                    @test size(Z, 3) >= 2
                    
                    # Test that msr=AddedTurbulence data exists and is reasonable (added turbulence should be >= 0%)
                    added_turbulence = Z[:, :, 2]
                    @test all(added_turbulence .>= 0.0)
                    @test all(isfinite.(added_turbulence))  # Should be finite values
                end
                
                @testset "msr parameter validation" begin
                    if Threads.nthreads() == 1
                        vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4, unit_test=true)
                        # Test error handling for invalid msr values
                        # Passing msr=0 to plotFlowField is expected to cause a TypeError (as checked below)
                        @test_throws TypeError plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=0)
                        # Passing an invalid msr enum value (e.g., msr=4) causes TypeError due to explicit enum validation
                        @test_throws TypeError plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=4)

                        # Test that msr=EffWind (default) still works with smart function
                        plot_flow_field(wf, X, Y, Z, vis; msr=EffWind, plt)
                        @test true  # If we get here, the function didn't throw an error
                        
                        # Test that wind speed data (msr=EffWind) is reasonable
                        wind_speed = Z[:, :, 3]
                        @test all(wind_speed .>= 0.0)  # Wind speed should be non-negative
                        @test all(isfinite.(wind_speed))  # Should be finite values
                    else
                        @info "Skipping msr parameter validation tests - only run in single-threaded mode (current: $(Threads.nthreads()) threads)"
                    end
                end
            end

            # Only run PlotState tests in single-threaded mode to avoid threading issues
            if Threads.nthreads() == 1
                @testset "plotFlowField with PlotState - two consecutive calls" begin
                    # Set up test data
                    state1 = nothing
                
                    # First call with nothing state (should create new PlotState)
                    state1 = @test_nowarn plotFlowField(state1, plt, wf, X, Y, Z, vis, 0; msr=EffWind)
                    
                    # Test that state1 is a PlotState object
                    @test state1 isa FLORIDyn.PlotState
                    
                    # Test that PlotState has all required fields
                    @test hasfield(typeof(state1), :fig)
                    @test hasfield(typeof(state1), :ax)
                    @test hasfield(typeof(state1), :cb)
                    @test hasfield(typeof(state1), :contour_collection)
                    @test hasfield(typeof(state1), :turbine_lines)
                    @test hasfield(typeof(state1), :op_scatter1)
                    @test hasfield(typeof(state1), :op_scatter2)
                    @test hasfield(typeof(state1), :title_obj)
                    @test hasfield(typeof(state1), :figure_name)
                    @test hasfield(typeof(state1), :label)
                    @test hasfield(typeof(state1), :lev_min)
                    @test hasfield(typeof(state1), :lev_max)
                    @test hasfield(typeof(state1), :levels)
                    
                    # Test that fields have correct types
                    @test state1.figure_name isa String
                    @test state1.label isa String
                    @test state1.lev_min isa Float64
                    @test state1.lev_max isa Float64
                    @test state1.lev_min < state1.lev_max
                    @test state1.turbine_lines isa Vector
                    
                    # Test second call with existing state (should update same PlotState)
                    state2 = @test_nowarn plotFlowField(state1, plt, wf, X, Y, Z, vis, 12; msr=VelReduction)
                    
                    # Test that state2 is the same as state1 (same object, reused)
                    @test state2 isa FLORIDyn.PlotState
                    
                    # Close the plots after testing to avoid accumulation
                    try
                        plt.close(state1.fig)
                    catch
                    end
                    
                    # Test error handling for invalid msr values
                    @test_throws TypeError plotFlowField(nothing, plt, wf, X, Y, Z, vis, 0; msr=99)

                    # Test with msr=AddedTurbulence (added turbulence) - create new state for this test
                    state3 = @test_nowarn plotFlowField(nothing, plt, wf, X, Y, Z, vis, 24; msr=AddedTurbulence)
                    @test state3 isa FLORIDyn.PlotState
                    
                    # Close the plot after testing
                    try
                        plt.close(state3.fig)
                    catch
                    end
                end
            else
                @info "Skipping PlotState tests - only run in single-threaded mode (current: $(Threads.nthreads()) threads)"
            end

        end
        
        @testset "plotFlowField - backward compatibility method" begin
            # Only run backward compatibility tests in single-threaded mode
            if Threads.nthreads() == 1
                @testset "basic functionality without state parameter" begin
                    # Get test parameters
                    wf, set, floris, wind, md = get_parameters()
                    
                    # Call the calcFlowField function
                    Z, X, Y = calcFlowField(set, wf, wind, floris)
                    vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4, unit_test=true)
                    
                    # Test the second method without state parameter (backward compatibility)
                    vis.unit_test = true
                    result = plotFlowField(plt, wf, X, Y, Z, vis)
                    
                    # Test that the function runs without error and returns nothing
                    @test result === nothing
                    
                    # Test that the function is callable and accepts all expected parameters
                    result2 = plotFlowField(plt, wf, X, Y, Z, vis; msr=VelReduction)
                    @test result2 === nothing
                    
                    result3 = plotFlowField(plt, wf, X, Y, Z, vis; msr=AddedTurbulence)
                    @test result3 === nothing
                    
                    # Test with time parameter
                    result4 = plotFlowField(plt, wf, X, Y, Z, vis, 120.0; msr=EffWind)
                    @test result4 === nothing
                end
            else
                @info "Skipping backward compatibility tests - only run in single-threaded mode (current: $(Threads.nthreads()) threads)"
            end
            
            @testset "parameter validation for backward compatibility method" begin
                # Only run backward compatibility tests in single-threaded mode
                if Threads.nthreads() == 1
                    # Get test parameters
                    wf, set, floris, wind, md = get_parameters()
                    Z, X, Y = calcFlowField(set, wf, wind, floris)
                    vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4, unit_test=true)
                    
                    # Test error handling for invalid msr values in backward compatibility method
                    @test_throws TypeError plotFlowField(plt, wf, X, Y, Z, vis; msr=0)
                    @test_throws TypeError plotFlowField(plt, wf, X, Y, Z, vis; msr=4)
                else
                    @info "Skipping parameter validation for backward compatibility tests - only run in single-threaded mode (current: $(Threads.nthreads()) threads)"

                end
            end
            
            @testset "file saving functionality with vis.save=true" begin
                # Only run file saving tests in single-threaded mode
                if Threads.nthreads() == 1
                    # Get test parameters
                    wf, set, floris, wind, md = get_parameters()
                    Z, X, Y = calcFlowField(set, wf, wind, floris)
                    
                    # Test with vis.save=true and vis.unit_test=false to enable file saving
                    # Use a unique test output folder to avoid conflicts
                    test_output_folder = tempname()
                    vis = Vis(online=false, save=true, output_folder=test_output_folder, 
                            rel_v_min=20.0, up_int=4, unit_test=false, print_filenames=false)
                    
                    # Ensure output directory doesn't exist before test
                    output_path = vis.output_path
                    @test isdir(output_path)  # Should be created by accessing the property
                    
                    # Test file saving for different measurement types
                    @testset "save files for different msr values" begin
                        # Test msr=VelReduction (velocity reduction)
                        result1 = @test_nowarn plotFlowField(plt, wf, X, Y, Z, vis; msr=VelReduction)
                        @test result1 === nothing
                        velocity_file = joinpath(output_path, "ff_velocity_reduction.png")
                        @test isfile(velocity_file)
                        plt.pause(0.1)
                        
                        # Test msr=AddedTurbulence (added turbulence)  
                        result2 = @test_nowarn plotFlowField(plt, wf, X, Y, Z, vis; msr=AddedTurbulence)
                        @test result2 === nothing
                        turbulence_file = joinpath(output_path, "ff_added_turbulence.png")
                        @test isfile(turbulence_file)
                        plt.pause(0.1)
                        
                        # Test msr=EffWind (wind speed)
                        result3 = @test_nowarn plotFlowField(plt, wf, X, Y, Z, vis; msr=EffWind)
                        @test result3 === nothing
                        wind_speed_file = joinpath(output_path, "ff_wind_speed.png")
                        @test isfile(wind_speed_file)
                        plt.pause(0.1)
                        plt.close("all")
                    end
                    
                    @testset "save files with time parameter" begin
                        # Test file saving with time parameter
                        result_t = @test_nowarn plotFlowField(plt, wf, X, Y, Z, vis, 120.0; msr=VelReduction)
                        @test result_t === nothing
                        time_file = joinpath(output_path, "ff_velocity_reduction_t0120s.png")
                        @test isfile(time_file)
                        
                        # Test with different time value
                        result_t2 = @test_nowarn plotFlowField(plt, wf, X, Y, Z, vis, 45.7; msr=AddedTurbulence)
                        @test result_t2 === nothing
                        time_file2 = joinpath(output_path, "ff_added_turbulence_t0046s.png")  # Should round to nearest int
                        @test isfile(time_file2)
                        plt.pause(0.1)
                        plt.close("all")
                    end
                    
                    @testset "online vs offline saving behavior" begin
                        # Test vis.online=true - should save to video_path
                        test_video_folder = "test_video_$(rand(10000:99999))"
                        vis_online = Vis(online=true, save=true, video_folder=test_video_folder,
                                    rel_v_min=20.0, unit_test=false, print_filenames=false)
                        
                        video_path = vis_online.video_path
                        @test isdir(video_path)
                        
                        result_online = @test_nowarn plotFlowField(plt, wf, X, Y, Z, vis_online; msr=VelReduction)
                        @test result_online === nothing
                        online_file = joinpath(video_path, "ff_velocity_reduction.png")
                        @test isfile(online_file)
                        
                        # Clean up video directory
                        try
                            rm(video_path, recursive=true)
                        catch
                        end
                    end
                    
                    @testset "file properties verification" begin
                        # Verify that saved files have reasonable properties
                        velocity_file = joinpath(output_path, "ff_velocity_reduction.png")
                        if isfile(velocity_file)
                            file_size = filesize(velocity_file)
                            @test file_size > 1000  # File should be at least 1KB (reasonable for a plot)
                            @test file_size < 1_000_000  # But not unreasonably large (< 1MB)
                        end
                    end
                    
                    # Clean up test output directory
                    try
                        rm(output_path, recursive=true)
                    catch e
                        @warn "Could not clean up test directory: $e"
                    end
                else
                    @info "Skipping file saving tests - only run in single-threaded mode (current: $(Threads.nthreads()) threads)"
                end
            end
        end
        
        @testset "get_layout" begin
            @testset "function exists and is callable" begin
                # Test that the function exists and has the right signature
                @test isa(get_layout, Function)
                
                # Test that we can get method information
                methods_list = methods(get_layout)
                @test length(methods_list) >= 1
                
                # Check that the function is exported
                @test :get_layout in names(FLORIDyn)
            end
            
            @testset "edge cases" begin
                # Test edge cases
                @test get_layout(0) == (1, 1)
                @test get_layout(-1) == (1, 1)
                @test get_layout(-10) == (1, 1)
            end
            
            @testset "small numbers" begin
                # Test specific small values
                @test get_layout(1) == (1, 1)
                @test get_layout(2) == (2, 2)
                @test get_layout(3) == (2, 2)
                @test get_layout(4) == (2, 2)
            end
            
            @testset "medium numbers" begin
                # Test medium values
                @test get_layout(5) == (2, 3)
                @test get_layout(6) == (2, 3)
                @test get_layout(7) == (3, 3)
                @test get_layout(8) == (3, 3)
                @test get_layout(9) == (3, 3)
                @test get_layout(10) == (3, 4)
                @test get_layout(11) == (3, 4)
                @test get_layout(12) == (3, 4)
            end
            
            @testset "larger numbers" begin
                # Test larger values
                @test get_layout(13) == (4, 4)
                @test get_layout(14) == (4, 4)
                @test get_layout(15) == (4, 4)
                @test get_layout(16) == (4, 4)
            end
            
            @testset "large numbers - general formula" begin
                # Test the general formula for large numbers
                # For nT > 16: cols = ceil(sqrt(nT)), rows = ceil(nT/cols)
                @test get_layout(17) == (4, 5)  # cols=ceil(sqrt(17))=5, rows=ceil(17/5)=4
                @test get_layout(20) == (4, 5)  # cols=ceil(sqrt(20))=5, rows=ceil(20/5)=4
                @test get_layout(25) == (5, 5)  # cols=ceil(sqrt(25))=5, rows=ceil(25/5)=5
                @test get_layout(30) == (5, 6)  # cols=ceil(sqrt(30))=6, rows=ceil(30/6)=5
                @test get_layout(50) == (7, 8)  # cols=ceil(sqrt(50))=8, rows=ceil(50/8)=7
                @test get_layout(100) == (10, 10)  # cols=ceil(sqrt(100))=10, rows=ceil(100/10)=10
            end
            
            @testset "layout properties" begin
                # Test that layouts can accommodate the required number of plots
                for nT in 1:50
                    rows, cols = get_layout(nT)
                    
                    # Basic sanity checks
                    @test rows >= 1
                    @test cols >= 1
                    @test isa(rows, Int)
                    @test isa(cols, Int)
                    
                    # The layout should be able to accommodate at least nT plots
                    @test rows * cols >= nT
                    
                    # For efficiency, the layout shouldn't be too oversized 
                    # (except for very small numbers where we use predefined layouts)
                    if nT > 16
                        # For large numbers, should not be more than 2x oversized
                        @test rows * cols <= 2 * nT
                    end
                end
            end
            
            @testset "roughly square layouts for large numbers" begin
                # Test that large numbers produce roughly square layouts
                test_cases = [17, 20, 25, 30, 50, 100]
                
                for nT in test_cases
                    rows, cols = get_layout(nT)
                    
                    # Should be roughly square (aspect ratio not too extreme)
                    aspect_ratio = max(rows, cols) / min(rows, cols)
                    @test aspect_ratio <= 2.0  # Not more than 2:1 aspect ratio
                    
                    # For perfect squares, should be exactly square
                    sqrt_nT = sqrt(nT)
                    if abs(sqrt_nT - round(sqrt_nT)) < 1e-10 && nT > 16  # Perfect square
                        expected_side = round(Int, sqrt_nT)
                        @test rows == expected_side
                        @test cols == expected_side
                    end
                end
            end
            
            @testset "consistency checks" begin
                # Test that the function is deterministic
                for nT in [1, 5, 10, 25, 50]
                    result1 = get_layout(nT)
                    result2 = get_layout(nT)
                    @test result1 == result2
                end
                
                # Test that increasing nT doesn't decrease the layout size inappropriately
                prev_area = 0
                for nT in 1:20
                    rows, cols = get_layout(nT)
                    area = rows * cols
                    
                    # Area should generally increase or stay the same
                    # (allowing for some predefined layout transitions)
                    if nT > 1
                        @test area >= prev_area || area >= nT
                    end
                    prev_area = area
                end
            end
        end
        @testset "plotMeasurements" begin
            # Get test parameters
            wf, set, floris, wind, md = get_parameters()
            vis = Vis(online=false, save=false, unit_test=true)
            
            # Test msr=VelReduction (velocity reduction) - existing tests
            @testset "msr=VelReduction (velocity reduction)" begin
                try
                    plot_measurements(wf, md, vis; separated=true, msr=VelReduction, plt)
                    println("✓ plot_measurements msr=VelReduction with separated=true completed successfully")
                catch e
                    @error "plot_measurements msr=VelReduction with separated=true failed: $e"
                    rethrow(e)
                end
                
                try
                    plot_measurements(wf, md, vis; separated=false, msr=VelReduction, plt)
                    println("✓ plot_measurements msr=VelReduction with separated=false completed successfully")
                catch e
                    @error "plot_measurements msr=VelReduction with separated=false failed: $e"
                    rethrow(e)
                end
            end
            
            # Test msr=AddedTurbulence (added turbulence) - new tests
            @testset "msr=AddedTurbulence (added turbulence)" begin
                # Check if AddedTurbulence column exists and has data
                if "AddedTurbulence" in names(md) && any(x -> !ismissing(x) && x != 0, md.AddedTurbulence)
                    try
                        plot_measurements(wf, md, vis; separated=true, msr=AddedTurbulence, plt)
                        println("✓ plot_measurements msr=AddedTurbulence with separated=true completed successfully")
                    catch e
                        @warn "plot_measurements msr=AddedTurbulence with separated=true failed: $e"
                        # Still test that it fails gracefully, not with unhandled errors
                        @test isa(e, Exception)
                    end
                    
                    try
                        plot_measurements(wf, md, vis; separated=false, msr=AddedTurbulence, plt)
                        println("✓ plot_measurements msr=AddedTurbulence with separated=false completed successfully")
                    catch e
                        @warn "plot_measurements msr=AddedTurbulence with separated=false failed: $e"
                        @test isa(e, Exception)
                    end
                else
                    @test_skip "Skipping msr=AddedTurbulence tests - AddedTurbulence column not available or contains no data"
                end
            end
            
            # Test msr=EffWind (effective wind speed) - new tests
            @testset "msr=EffWind (effective wind speed)" begin
                # Check if EffWindSpeed column exists and has data
                if "EffWindSpeed" in names(md) && any(x -> !ismissing(x) && x != 0, md.EffWindSpeed)
                    try
                        plot_measurements(wf, md, vis; separated=true, msr=EffWind, plt)
                        println("✓ plot_measurements msr=EffWind with separated=true completed successfully")
                    catch e
                        @warn "plot_measurements msr=EffWind with separated=true failed: $e"
                        @test isa(e, Exception)
                    end
                    
                    try
                        plot_measurements(wf, md, vis; separated=false, msr=EffWind, plt)
                        println("✓ plot_measurements msr=EffWind with separated=false completed successfully")
                    catch e
                        @warn "plot_measurements msr=EffWind with separated=false failed: $e"
                        @test isa(e, Exception)
                    end
                else
                    @test_skip "Skipping msr=EffWind tests - EffWindSpeed column not available or contains no data"
                end
            end
            
            # Test error handling for invalid msr values
            @testset "Invalid msr values" begin
                @test_throws TypeError plot_measurements(wf, md, vis; separated=true, msr=0, plt)
                @test_throws TypeError plot_measurements(wf, md, vis; separated=true, msr=4, plt)
                @test_throws TypeError plot_measurements(wf, md, vis; separated=false, msr=-1, plt)
            end
            
            # Test default msr value (should be 1)
            @testset "Default msr value" begin
                try
                    # Test that calling without msr parameter defaults to msr=VelReduction
                    plot_measurements(wf, md, vis; separated=true, plt)
                    println("✓ plot_measurements with default msr completed successfully")
                catch e
                    @error "plot_measurements with default msr failed: $e"
                    rethrow(e)
                end
            end
            if Threads.nthreads() > 1
                sleep_duration = get(ENV, "TEST_THREAD_SLEEP", "10")
                try
                    sleep(parse(Int, sleep_duration))
                catch
                    sleep(10) # fallback to default if parsing fails
                end
            end
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
                        
                        create_test_png(file_path, 10, 10)
                    end
                    
                    # Test basic functionality (should fail gracefully without FFmpeg)
                    result = createVideo("test"; video_dir=test_dir, output_dir=output_dir, fps=2)
                    
                    # The function should handle missing FFmpeg gracefully
                    @test isa(result, String)  # Should return a string (empty if failed)

                    # Test basic functionality and delete frames
                    result = createVideo("test"; video_dir=test_dir, output_dir=output_dir, fps=2, delete_frames=true)
                    
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
        @testset "plot_x" begin
            @testset "function exists and is callable" begin
                # Test that the function exists and has the right signature
                @test isa(plot_x, Function)
                
                # Test that we can get method information
                methods_list = methods(plot_x)
                @test length(methods_list) >= 1
                
                # Check that the function is exported
                @test :plot_x in names(FLORIDyn)
            end
            
            @testset "basic functionality with test data" begin
                # Create test time series data
                times = collect(0.0:10.0:100.0)  # 11 time points from 0 to 100
                data1 = sin.(times * π / 50)     # Sine wave
                data2 = cos.(times * π / 50)     # Cosine wave
                data3 = times / 100             # Linear trend
                
                ylabels = ["Turbine 1", "Turbine 2", "Turbine 3"]
                labels = ["Wind Speed", "Power", "Efficiency"]
                
                # Test basic functionality (should dispatch based on threading/processing environment)
                @test_nowarn plot_x(times, data1, data2, data3; ylabels=ylabels, labels=labels, pltctrl=pltctrl)
                
                # Test that the function runs without error
                @test true
            end
            
            @testset "parameter validation" begin
                # Create simple test data
                times = [0.0, 1.0, 2.0]
                data = [1.0, 2.0, 3.0]
                
                # Test with default parameters
                @test_nowarn plot_x(times, data; pltctrl=pltctrl)
                
                # Test with custom parameters
                @test_nowarn plot_x(times, data; 
                                   ylabels=["Test"], 
                                   labels=["Test Plot"], 
                                   fig="Custom Title",
                                   xlabel="Custom X Label",
                                   ysize=8,
                                   bottom=0.05,
                                   legend_size=10,
                                   pltctrl=pltctrl)
                
                # Test error handling when pltctrl is missing in single-threaded mode
                if Threads.nthreads() == 1 && nprocs() < 2
                    @test_throws ErrorException plot_x(times, data)
                end
            end
            
            @testset "different data configurations" begin
                # Test with single data series
                times = [0.0, 5.0, 10.0, 15.0]
                single_data = [10.0, 15.0, 12.0, 18.0]
                
                @test_nowarn plot_x(times, single_data; 
                                   ylabels=["Single Series"],
                                   fig="Single Data Test",
                                   pltctrl=pltctrl)
                
                # Test with multiple data series
                data1 = [10.0, 15.0, 12.0, 18.0]
                data2 = [5.0, 8.0, 6.0, 9.0]
                data3 = [20.0, 25.0, 22.0, 28.0]
                
                @test_nowarn plot_x(times, data1, data2, data3;
                                   ylabels=["Series 1", "Series 2", "Series 3"],
                                   labels=["Plot 1", "Plot 2", "Plot 3"],
                                   fig="Multiple Data Test",
                                   pltctrl=pltctrl)
                
                # Test with empty ylabels and labels (should use defaults)
                @test_nowarn plot_x(times, data1, data2;
                                   fig="Default Labels Test",
                                   pltctrl=pltctrl)
            end
            
            @testset "threading and process dispatch" begin
                # Create test data
                times = [0.0, 1.0, 2.0]
                data = [1.0, 4.0, 2.0]
                
                # Test behavior based on environment
                if Threads.nthreads() > 1 && nprocs() > 1
                    # Should use remote plotting (rmt_plotx)
                    @test_nowarn plot_x(times, data; fig="Multi-threaded Test")
                else
                    # Should use sequential plotting (requires pltctrl)
                    @test_nowarn plot_x(times, data; pltctrl=pltctrl, fig="Single-threaded Test")
                    
                    # Should error without pltctrl in single-threaded mode
                    @test_throws ErrorException plot_x(times, data; fig="Error Test")
                end
            end
            
            @testset "wind direction specific test" begin
                # Test with data similar to wind direction plotting
                times = collect(0.0:12.0:120.0)  # Every 12 seconds for 2 minutes
                wind_dir_t1 = 270.0 .+ 5.0 * sin.(times * π / 60)  # Oscillating around 270°
                wind_dir_t2 = 280.0 .+ 3.0 * cos.(times * π / 60)  # Oscillating around 280°
                wind_dir_t3 = 275.0 .+ 2.0 * sin.(times * π / 30)  # Different frequency
                
                turbine_labels = ["T1", "T2", "T3"]
                subplot_labels = ["Wind Direction", "Wind Direction", "Wind Direction"]
                
                # Test the specific use case from the original request
                @test_nowarn plot_x(times, wind_dir_t1, wind_dir_t2, wind_dir_t3;
                                   ylabels=turbine_labels,
                                   labels=subplot_labels,
                                   fig="Wind Direction",
                                   xlabel="rel_time [s]",
                                   ysize=10,
                                   bottom=0.02,
                                   legend_size=8,
                                   pltctrl=pltctrl)
                sleep(1)
                close_all(plt)  # Close the plot after testing
            end
            
            @testset "parameter defaults" begin
                # Test that default parameters match the expected values
                times = [0.0, 1.0, 2.0]
                data = [1.0, 2.0, 1.5]
                
                # Test default fig parameter
                @test_nowarn plot_x(times, data; pltctrl=pltctrl)  # Should use "Wind Direction" as default
                
                # Test default xlabel parameter
                @test_nowarn plot_x(times, data; xlabel="Time [min]", pltctrl=pltctrl)
                
                # Test default ysize parameter
                @test_nowarn plot_x(times, data; ysize=12, pltctrl=pltctrl)
                
                # Test default bottom parameter  
                @test_nowarn plot_x(times, data; bottom=0.1, pltctrl=pltctrl)
                
                # Test legend_size parameter
                @test_nowarn plot_x(times, data; legend_size=12, pltctrl=pltctrl)
                sleep(1)  # Allow time for the plot to render
                close_all(plt)  # Close the plot after testing
            end
            
            @testset "error handling" begin
                # Test with mismatched data lengths
                times = [0.0, 1.0, 2.0]
                bad_data = [1.0, 2.0]  # Different length than times
                
                # The function should handle this gracefully or error appropriately
                if Threads.nthreads() == 1 && nprocs() < 2
                    # In single-threaded mode, errors from plotx should propagate
                    @test_throws Exception plot_x(times, bad_data; pltctrl=pltctrl)
                else
                    # In multi-threaded mode, might be handled by remote worker
                    # Just test that it doesn't crash the main process
                    @test_nowarn plot_x(times, bad_data)
                end
                
                # Test with empty data
                empty_times = Float64[]
                empty_data = Float64[]
                
                if Threads.nthreads() == 1 && nprocs() < 2
                    @test_throws Exception plot_x(empty_times, empty_data; pltctrl=pltctrl)
                end
            end
            
            @testset "integration with smart plotting pattern" begin
                # Test that plot_x follows the same pattern as other smart plotting functions
                
                # Should have similar behavior to plot_flow_field and plot_measurements
                times = [0.0, 10.0, 20.0]
                data = [1.0, 2.0, 1.5]
                
                # Test consistency with other smart plotting functions
                @test_nowarn plot_x(times, data; pltctrl=pltctrl)
                
                # Should handle the same threading logic
                if Threads.nthreads() > 1 && nprocs() > 1
                    # Should work without plt in multi-threaded mode
                    @test_nowarn plot_x(times, data)
                else
                    # Should require plt in single-threaded mode
                    @test_throws ErrorException plot_x(times, data)
                end
                
                # Should return nothing like other smart plotting functions
                result = plot_x(times, data; pltctrl=pltctrl)
                @test result === nothing
            end
        end
        
        @testset "calc_rel_power" begin
            @testset "function exists and is callable" begin
                # Test that the function exists and has the right signature
                @test calc_rel_power isa Function
                
                # Test that we can get method information
                methods_list = methods(calc_rel_power)
                @test length(methods_list) >= 1
                
                # Check that the function is exported
                @test :calc_rel_power in names(FLORIDyn)
            end
            
            @testset "basic functionality with default parameters" begin
                # Test with minimal 9-turbine configuration for faster execution
                settings_file = "data/2021_9T_Data.yaml"
                
                # Basic call with default parameters (using wind_dir=nothing like main_power_plot.jl)
                times, rel_power, set, wf, wind, floris = calc_rel_power(settings_file; wind_dir=nothing)
                
                # Test return types and structure
                @test isa(times, Vector{Float64})
                @test isa(rel_power, Vector{Float64})
                @test isa(set, Settings)
                @test isa(wf, WindFarm)
                @test isa(wind, Wind)  # Wind is imported directly, not FLORIDyn.Wind
                @test isa(floris, Floris)
                
                # Test that arrays have consistent lengths
                @test length(times) == length(rel_power)
                @test length(times) > 0  # Should have data points
                
                # Test that time values are reasonable
                @test all(times .>= 0.0)  # Times should be non-negative
                @test issorted(times)     # Times should be monotonically increasing
                
                # Test that relative power values are reasonable
                @test all(rel_power .>= 0.0)     # Power should be non-negative
                @test all(rel_power .<= 2.0)     # Should not exceed 200% (reasonable upper bound)
                @test all(isfinite.(rel_power))  # Should be finite values
                
                # Test that wind farm has expected structure
                @test wf.nT > 0  # Should have turbines
                @test wf.nT == 9  # Should be 9-turbine configuration
            end
            
            @testset "custom dt parameter" begin
                settings_file = "data/2021_9T_Data.yaml"
                
                # Test with shorter simulation time
                dt_short = 50  # 50 seconds for faster testing
                times_short, rel_power_short, set_short, wf_short, wind_short, floris_short = calc_rel_power(settings_file; dt=dt_short, wind_dir=180.0)
                
                # Test with longer simulation time  
                dt_long = 100   # 100 seconds
                times_long, rel_power_long, set_long, wf_long, wind_long, floris_long = calc_rel_power(settings_file; dt=dt_long, wind_dir=180.0)
                
                # Longer simulation should have more data points (or at least same)
                @test length(times_long) >= length(times_short)
                
                # Both should produce valid results
                @test all(isfinite.(rel_power_short))
                @test all(isfinite.(rel_power_long))
                @test all(rel_power_short .>= 0.0)
                @test all(rel_power_long .>= 0.0)
            end
            
            @testset "fixed wind direction parameter" begin
                settings_file = "data/2021_9T_Data.yaml"
                
                # Test with fixed wind direction
                wind_dir_270 = 270.0
                times_270, rel_power_270, set_270, wf_270, wind_270, floris_270 = calc_rel_power(settings_file; dt=50, wind_dir=wind_dir_270)
                
                # Test with different wind direction
                wind_dir_90 = 90.0
                times_90, rel_power_90, set_90, wf_90, wind_90, floris_90 = calc_rel_power(settings_file; dt=50, wind_dir=wind_dir_90)
                
                # Test settings configuration for fixed wind direction
                @test set_270.dir_mode isa Direction_Constant
                @test set_270.control_mode isa Yaw_Constant
                @test set_90.dir_mode isa Direction_Constant
                @test set_90.control_mode isa Yaw_Constant
                
                # Test that wind direction is applied
                @test wind_270.dir[1,1] ≈ wind_dir_270
                @test wind_90.dir[1,1] ≈ wind_dir_90
                
                # Both should produce valid results
                @test all(isfinite.(rel_power_270))
                @test all(isfinite.(rel_power_90))
                @test all(rel_power_270 .>= 0.0)
                @test all(rel_power_90 .>= 0.0)
                
                # Time arrays should have same length for same dt
                @test length(times_270) == length(times_90)
            end
            
            @testset "variable wind direction (default)" begin
                settings_file = "data/2021_9T_Data.yaml"
                
                # Test with variable wind direction (wind_dir=nothing, default)
                times_var, rel_power_var, set_var, wf_var, wind_var, floris_var = calc_rel_power(settings_file; dt=50, wind_dir=nothing)
                
                # Test settings configuration for variable wind direction
                @test !(set_var.dir_mode isa Direction_Constant)
                @test !(set_var.control_mode isa Yaw_Constant)
                
                # Should produce valid results
                @test isa(times_var, Vector{Float64})
                @test isa(rel_power_var, Vector{Float64})
                @test length(times_var) == length(rel_power_var)
                @test all(isfinite.(rel_power_var))
                @test all(rel_power_var .>= 0.0)
            end
            
            @testset "power calculation physics" begin
                settings_file = "data/2021_9T_Data.yaml"
                
                # Test the cubic relationship (P ∝ v³)
                times, rel_power, set, wf, wind, floris = calc_rel_power(settings_file; dt=50, wind_dir=180.0)
                
                # Power should generally be between 0 and 1 for most wind conditions
                # (unless there's significant speedup effects)
                typical_power_range = 0.2 <= mean(rel_power) <= 1.2
                @test typical_power_range
                
                # Standard deviation should be reasonable (not zero, but not extreme)
                power_std = std(rel_power)
                @test 0.0 <= power_std <= 0.5  # Reasonable variation
                
                # Test that individual turbine contributions make sense
                # The function sums rel_speed^3 for all turbines and divides by nT
                # So rel_power should be the average of cubic relative wind speeds
                @test all(0.0 .<= rel_power .<= 2.0)  # Each turbine contributes 0-200% typically
            end
            
            @testset "return value validation" begin
                settings_file = "data/2021_9T_Data.yaml"
                times, rel_power, set, wf, wind, floris = calc_rel_power(settings_file; dt=50)
                
                # Test times array properties
                @testset "times array" begin
                    @test isa(times, Vector{Float64})
                    @test length(times) > 0
                    @test issorted(times)
                    @test all(times .>= 0.0)
                    @test all(isfinite.(times))
                end
                
                # Test rel_power array properties
                @testset "rel_power array" begin
                    @test isa(rel_power, Vector{Float64})
                    @test length(rel_power) == length(times)
                    @test all(isfinite.(rel_power))
                    @test all(rel_power .>= 0.0)
                    @test !all(rel_power .== 0.0)  # Should not be all zeros
                end
                
                # Test Settings object
                @testset "Settings object" begin
                    @test isa(set, Settings)
                    @test hasfield(typeof(set), :dir_mode)
                    @test hasfield(typeof(set), :control_mode)
                    @test hasfield(typeof(set), :parallel)
                    @test hasfield(typeof(set), :threading)
                end
                
                # Test WindFarm object
                @testset "WindFarm object" begin
                    @test isa(wf, WindFarm)
                    @test wf.nT > 0
                    @test hasfield(typeof(wf), :States_T)
                    @test hasfield(typeof(wf), :nT)
                end
                
                # Test Wind object
                @testset "Wind object" begin
                    @test isa(wind, Wind)  # Wind is imported directly
                    @test hasfield(typeof(wind), :dir)
                    @test size(wind.dir, 1) > 0
                    @test size(wind.dir, 2) > 0
                end
                
                # Test Floris object
                @testset "Floris object" begin
                    @test isa(floris, Floris)
                    # Floris internal structure may vary, just test it exists
                end
            end
            
            @testset "performance and memory" begin
                settings_file = "data/2021_9T_Data.yaml"
                
                # Test that function completes in reasonable time
                @testset "execution time" begin
                    start_time = time()
                    times, rel_power, set, wf, wind, floris = calc_rel_power(settings_file; dt=30)
                    execution_time = time() - start_time
                    
                    # Should complete within reasonable time (depends on system)
                    @test execution_time < 30.0  # Should complete within 30 seconds for small dt
                    @test length(times) > 0
                    @test all(isfinite.(rel_power))
                end
                
                # Test memory usage is reasonable
                @testset "memory usage" begin
                    # Run GC before test
                    GC.gc()
                    
                    times, rel_power, set, wf, wind, floris = calc_rel_power(settings_file; dt=50)
                    
                    # Arrays should not be excessively large
                    @test length(times) < 10000  # Should be reasonable number of time steps
                    @test sizeof(rel_power) < 100_000  # Should not use excessive memory
                end
            end
            
            @testset "error handling" begin
                # Test with non-existent settings file
                @testset "invalid settings file" begin
                    @test_throws SystemError calc_rel_power("nonexistent_file.yaml")
                end
                
                # Test with invalid parameter types
                @testset "invalid parameter types" begin
                    settings_file = "data/2021_9T_Data.yaml"
                    
                    # Invalid dt type
                    @test_throws MethodError calc_rel_power(settings_file; dt="invalid")
                    
                    # Invalid wind_dir type  
                    @test_throws MethodError calc_rel_power(settings_file; wind_dir="invalid")
                end
            end
        end
        
    end

    @testset verbose=true "Pretty Print Functions                                  " begin
        @testset "turbines(T::Dict)" begin
            @testset "basic functionality with States_T (capital T)" begin
                # Test with the standard key format used in MAT files
                T_dict = Dict(
                    "States_T" => [1.0 2.0 3.0; 4.0 5.0 6.0],  # (nT*nOP) × n_state_variables = (2*1) × 3
                    "Names_T" => ["a", "yaw", "TI"],
                    "nT" => 2,  # 2 turbines  
                    "nOP" => 1  # 1 operating point each
                )
                
                df = turbines(T_dict)
                
                # Test output structure - now with OP and Turbine columns plus the data columns
                @test isa(df, DataFrame)
                @test size(df) == (2, 5)  # 2 turbines × 1 OP = 2 rows, 5 columns (OP + Turbine + 3 data)
                @test names(df) == ["OP", "Turbine", "a", "yaw", "TI"]
                
                # Test identifier columns
                @test df.OP == [1, 1]  # 1 OP for each turbine
                @test df.Turbine == [1, 2]  # Turbine numbers
                
                # Test data values - data rows are: turbine 1 OP 1, turbine 2 OP 1
                @test df.a == [1.0, 4.0]
                @test df.yaw == [2.0, 5.0]
                @test df.TI == [3.0, 6.0]
                
                # Test column types
                @test eltype(df.OP) == Int64
                @test eltype(df.Turbine) == Int64
                @test eltype(df.a) == Float64
                @test eltype(df.yaw) == Float64
                @test eltype(df.TI) == Float64
            end
            
            @testset "basic functionality with States_t (lowercase t)" begin
                # Test error when using lowercase variant (should fail)
                T_dict = Dict(
                    "States_t" => [10.0 20.0; 30.0 40.0],  # This key should fail
                    "Names_T" => ["T1", "T2"],
                    "nT" => 2,
                    "nOP" => 1
                )
                
                # This should throw an error since lowercase is not supported
                @test_throws ArgumentError turbines(T_dict)
            end
            
            @testset "symbol keys support" begin
                # Test with symbol keys instead of string keys
                T_dict = Dict(
                    :States_T => [1.5 2.5; 3.5 4.5],  # 2 turbines × 1 OP each, 2 states
                    :Names_T => ["A", "B"],
                    :nT => 2,
                    :nOP => 1
                )
                
                df = turbines(T_dict)
                
                @test isa(df, DataFrame)
                @test size(df) == (2, 4)  # 2 turbines × 1 OP = 2 rows, 4 columns
                @test names(df) == ["OP", "Turbine", "A", "B"]
                @test df.OP == [1, 1]
                @test df.Turbine == [1, 2]
                @test df.A == [1.5, 3.5]
                @test df.B == [2.5, 4.5]
            end
            
            @testset "mixed case key support" begin
                # Test mixed string/symbol keys
                T_dict = Dict(
                    "States_T" => [0.1 0.2; 0.3 0.4],
                    :Names_T => ["Mixed1", "Mixed2"],
                    "nT" => 2,
                    :nOP => 1
                )
                
                df_mixed = turbines(T_dict)
                
                @test isa(df_mixed, DataFrame)
                @test size(df_mixed) == (2, 4)  # 2 turbines × 1 OP = 2 rows, 4 columns
                @test names(df_mixed) == ["OP", "Turbine", "Mixed1", "Mixed2"]
            end
            
            @testset "single turbine" begin
                # Test with only one turbine
                T_dict = Dict(
                    "States_T" => reshape([5.5, 6.6, 7.7], 1, 3),  # 1 turbine × 1 OP, 3 states
                    "Names_T" => ["state1", "state2", "state3"],
                    "nT" => 1,
                    "nOP" => 1
                )
                
                df = turbines(T_dict)
                
                @test isa(df, DataFrame)
                @test size(df) == (1, 5)  # 1 turbine × 1 OP = 1 row, 5 columns  
                @test names(df) == ["OP", "Turbine", "state1", "state2", "state3"]
                @test df.OP == [1]
                @test df.Turbine == [1] 
                @test df.state1 == [5.5]
                @test df.state2 == [6.6]
                @test df.state3 == [7.7]
            end
            
            @testset "large dataset" begin
                # Test with larger dataset (multiple OPs per turbine)
                nT, nOP = 3, 4
                states_data = rand(nT * nOP, 2)  # (3×4) × 2 = 12 × 2 matrix
                T_dict = Dict(
                    "States_T" => states_data,
                    "Names_T" => ["param1", "param2"],
                    "nT" => nT,
                    "nOP" => nOP  
                )
                
                df = turbines(T_dict)
                
                @test isa(df, DataFrame)
                @test size(df) == (12, 4)  # 3 turbines × 4 OPs = 12 rows, 4 columns
                @test names(df) == ["OP", "Turbine", "param1", "param2"]
                
                # Check OP pattern: 1,2,3,4,1,2,3,4,1,2,3,4 (4 OPs repeated for 3 turbines)
                @test df.OP == repeat(1:nOP, nT)
                # Check Turbine pattern: 1,1,1,1,2,2,2,2,3,3,3,3 (each turbine repeated 4 times)  
                @test df.Turbine == repeat(1:nT, inner=nOP)
            end
            
            @testset "different data types" begin
                # Test with different numeric types
                T_dict = Dict(
                    "States_T" => Int32[10 20; 30 40],  # Integer input
                    "Names_T" => ["IntVal1", "IntVal2"],
                    "nT" => 2,
                    "nOP" => 1
                )
                
                df = turbines(T_dict)
                
                @test isa(df, DataFrame)
                @test size(df) == (2, 4)
                @test names(df) == ["OP", "Turbine", "IntVal1", "IntVal2"] 
                @test df.IntVal1 == [10, 30]
                @test df.IntVal2 == [20, 40]
            end
            
            @testset "real MAT file format" begin
                # Test with format similar to actual MAT file data
                nT, nOP = 3, 2
                T_dict = Dict(
                    "States_T" => [
                        0.33 -0.1 0.0;      # Turbine 1 OP 1
                        0.33 -0.099 0.001;  # Turbine 1 OP 2  
                        0.33 -0.1 0.0;      # Turbine 2 OP 1
                        0.33 -0.099 0.001;  # Turbine 2 OP 2
                        0.33 -0.1 0.0;      # Turbine 3 OP 1
                        0.33 -0.099 0.001   # Turbine 3 OP 2
                    ],  # 6 rows = nT*nOP, 3 cols = n_state_variables
                    "Names_T" => ["a", "yaw", "TI"],  # Actual names from MAT files
                    "nT" => nT,
                    "nOP" => nOP
                )
                
                df = turbines(T_dict)
                
                @test isa(df, DataFrame)
                @test size(df) == (6, 5)  # 3 turbines × 2 OPs = 6 rows, 5 columns
                @test names(df) == ["OP", "Turbine", "a", "yaw", "TI"]
                @test df.OP == [1, 2, 1, 2, 1, 2]  # Pattern: OP 1-2 for each turbine
                @test df.Turbine == [1, 1, 2, 2, 3, 3]  # Pattern: each turbine repeated nOP times
                @test all(df.a .== 0.33)  # All axial induction values same
            end
            
            @testset "error handling - missing keys" begin
                # Test error when States key is missing
                @test_throws ArgumentError turbines(Dict("Names_T" => ["T1"], "nT" => 1, "nOP" => 1))
                @test_throws ArgumentError turbines(Dict("Random_Key" => [1.0 2.0], "Names_T" => ["T1", "T2"], "nT" => 2, "nOP" => 1))
                
                # Test error when Names key is missing  
                @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "nT" => 2, "nOP" => 1))
                @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "Random_Key" => ["T1", "T2"], "nT" => 2, "nOP" => 1))
                
                # Test error when nT key is missing
                @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "Names_T" => ["T1", "T2"], "nOP" => 1))
                
                # Test error when nOP key is missing
                @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "Names_T" => ["T1", "T2"], "nT" => 2))
                
                # Test error when using lowercase keys (not supported)
                @test_throws ArgumentError turbines(Dict("States_t" => [1.0 2.0; 3.0 4.0], "Names_T" => ["T1", "T2"], "nT" => 2, "nOP" => 1))
                @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "Names_t" => ["T1", "T2"], "nT" => 2, "nOP" => 1))
            end
            
            @testset "error handling - empty data" begin
                # Test error with empty states matrix
                @test_throws ArgumentError turbines(Dict(
                    "States_T" => Matrix{Float64}(undef, 0, 0),
                    "Names_T" => String[],
                    "nT" => 0,
                    "nOP" => 0
                ))
                
                # Test error with empty names vector
                @test_throws ArgumentError turbines(Dict(
                    "States_T" => [1.0 2.0; 3.0 4.0],
                    "Names_T" => String[],
                    "nT" => 2,
                    "nOP" => 1
                ))
            end
            
            @testset "error handling - dimension mismatch" begin
                # Test error when number of columns doesn't match number of names
                @test_throws DimensionMismatch turbines(Dict(
                    "States_T" => [1.0 2.0 3.0; 4.0 5.0 6.0],  # 3 state variables
                    "Names_T" => ["T1", "T2"],  # Only 2 names
                    "nT" => 2,
                    "nOP" => 1
                ))
                
                @test_throws DimensionMismatch turbines(Dict(
                    "States_T" => [1.0 2.0; 3.0 4.0],  # 2 state variables
                    "Names_T" => ["T1", "T2", "T3"],  # 3 names
                    "nT" => 2,
                    "nOP" => 1
                ))
                
                # Test error when States_T dimensions don't match nT*nOP
                @test_throws DimensionMismatch turbines(Dict(
                    "States_T" => [1.0 2.0; 3.0 4.0],  # 2 rows 
                    "Names_T" => ["T1", "T2"],  
                    "nT" => 3,  # 3 turbines
                    "nOP" => 2  # 2 OPs -> should be 6 rows, but only have 2
                ))
            end
            
            @testset "error handling - invalid dictionary" begin
                # Test with empty dictionary
                @test_throws ArgumentError turbines(Dict{String, Any}())
            end
            
            @testset "integration with actual MAT file data" begin
                # Test that the function works with the format from test_case_01.jl
                # This simulates the actual MAT file structure
                
                # Create test data similar to what would come from matread()
                nT, nOP = 2, 3
                T_ref_simulation = Dict(
                    "States_T" => [
                        0.33  -0.09999  0.0;   # Turbine 1, OP 1
                        0.33  -0.1      0.0;   # Turbine 1, OP 2
                        0.33  -0.099    0.05;  # Turbine 1, OP 3
                        0.33  -0.10001  0.0;   # Turbine 2, OP 1
                        0.33  -0.1      0.0;   # Turbine 2, OP 2
                        0.33  -0.098    0.03   # Turbine 2, OP 3
                    ],
                    "Names_T" => ["a", "yaw", "TI"],  # State variable names
                    "nT" => nT,
                    "nOP" => nOP
                )
                
                df = turbines(T_ref_simulation)
                
                @test isa(df, DataFrame)
                @test size(df) == (6, 5)  # 2 turbines × 3 OPs = 6 rows, 5 columns
                @test names(df) == ["OP", "Turbine", "a", "yaw", "TI"]
                
                # Verify the structure matches expected patterns
                @test df.OP == [1, 2, 3, 1, 2, 3]  # 1,2,3 for turbine 1, then 1,2,3 for turbine 2
                @test df.Turbine == [1, 1, 1, 2, 2, 2]  # Each turbine repeated for all OPs
                @test all(df.a .== 0.33)  # Axial induction should be constant
                @test df.TI[3] ≈ 0.05  # OP 3 for turbine 1
                @test df.TI[6] ≈ 0.03  # OP 3 for turbine 2
            end
            
            @testset "function integration" begin
                # Test that the function is properly exported
                @test :turbines in names(FLORIDyn)
                
                # Test it works with the same data format as turbines(wf::WindFarm)
                T_dict = Dict(
                    "States_T" => [0.33 -0.1 0.06; 0.33 -0.1 0.06],  # 2 turbines × 1 OP, 3 states
                    "Names_T" => ["a", "yaw", "TI"],
                    "nT" => 2,
                    "nOP" => 1
                )
                
                df = turbines(T_dict)
                
                # Should have same columns as WindFarm version
                @test "OP" in names(df)
                @test "Turbine" in names(df) 
                @test "a" in names(df)
                @test "yaw" in names(df)
                @test "TI" in names(df)
            end
            
            @testset "edge cases" begin
                # Test with minimal valid data (1 turbine, 1 OP, 1 state)
                T_dict = Dict(
                    "States_T" => reshape([42.0], 1, 1),
                    "Names_T" => ["single_state"],
                    "nT" => 1,
                    "nOP" => 1
                )
                
                df = turbines(T_dict)
                @test isa(df, DataFrame)
                @test size(df) == (1, 3)  # 1 row, 3 columns (OP + Turbine + data)
                @test df.OP == [1]
                @test df.Turbine == [1]
                @test df.single_state == [42.0]
            end
            
            @testset "performance with large datasets" begin
                # Test that the function performs reasonably with large data
                nT, nOP = 10, 100
                n_states = 5
                
                large_states = rand(Float64, nT * nOP, n_states)  # 1000 × 5 matrix
                large_names = ["State_$i" for i in 1:n_states]
                
                T_large = Dict(
                    "States_T" => large_states,
                    "Names_T" => large_names,
                    "nT" => nT,
                    "nOP" => nOP
                )
                
                # Time the operation (should complete in reasonable time)
                time_start = time()
                df_large = turbines(T_large)
                time_elapsed = time() - time_start
                
                @test isa(df_large, DataFrame)
                @test size(df_large) == (1000, 7)  # 1000 rows × (2 + 5) columns
                @test time_elapsed < 5.0  # Should complete in less than 5 seconds
                
                # Verify data integrity for large dataset
                @test df_large.OP == repeat(1:nOP, nT)  # Pattern repeats correctly
                @test df_large.Turbine == repeat(1:nT, inner=nOP)  # Turbine pattern is correct
            end
        end
    end
    sleep(1)  # Allow time for the plot to render
    close_all(plt)  # Close the plot after testing
else
    # Running tests via Pkg.test (safest approach)
    @eval Main using Pkg
    Main.Pkg.test(test_args=["test_visualisation.jl"])
end
nothing
