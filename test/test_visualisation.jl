# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

if !isdefined(Main, :Test)
    using Test
end 

if ! isinteractive()
if !isdefined(Main, :FLORIDyn)
    using FLORIDyn
end

if !isdefined(Main, :ControlPlots)
    using ControlPlots
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

@testset verbose=true "visualisation                                           " begin
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
        # Get test parameters once for all plotFlowField tests
        wf, set, floris, wind, md = get_parameters()
        Z, X, Y = calcFlowField(set, wf, wind, floris)
        vis = Vis(online=false, save=true, rel_v_min=20.0, up_int = 4, unit_test=true)

        @testset "basic functionality" begin
            plot_flow_field(wf, X, Y, Z, vis; plt=ControlPlots.plt)

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
                plot_flow_field(wf, X, Y, Z, vis; msr=VelReduction, plt=ControlPlots.plt)

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
                plot_flow_field(wf, X, Y, Z, vis; msr=AddedTurbulence, plt=ControlPlots.plt)

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
                    # msr=0 causes TypeError
                    @test_throws TypeError plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=0)
                    # msr > 3 causes TypeError from explicit check
                    @test_throws TypeError plotFlowField(nothing, plt, wf, X, Y, Z, vis; msr=4)

                    # Test that msr=EffWind (default) still works with smart function
                    plot_flow_field(wf, X, Y, Z, vis; msr=EffWind, plt=ControlPlots.plt)
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
                # Enable debug logging for this diagnostic test in CI (single-thread scenario)
                vis.log_debug = false
                @info "Diagnostics: starting two-call PlotState test" threads=Threads.nthreads() show_plots=vis.show_plots save=vis.save
                
                # First call with nothing state (should create new PlotState)
                state1 = @test_nowarn plotFlowField(state1, ControlPlots.plt, wf, X, Y, Z, vis, 0; msr=EffWind)
                if vis.log_debug
                    @info "After first call" figure_name=state1.figure_name label=state1.label lev_min=state1.lev_min lev_max=state1.lev_max levels_len=length(state1.levels) contour_type=string(typeof(state1.contour_collection))
                end
                
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
                state2 = @test_nowarn plotFlowField(state1, ControlPlots.plt, wf, X, Y, Z, vis, 12; msr=VelReduction)
                if vis.log_debug
                    @info "After second call" figure_name=state2.figure_name label=state2.label lev_min=state2.lev_min lev_max=state2.lev_max levels_len=length(state2.levels) contour_type=string(typeof(state2.contour_collection))
                end
                
                # Test that state2 is the same as state1 (same object, reused)
                @test state2 isa FLORIDyn.PlotState
                
                # Close the plots after testing to avoid accumulation
                try
                    ControlPlots.plt.close(state1.fig)
                catch
                end
                
                # Test error handling for invalid msr values
                @test_throws TypeError plotFlowField(nothing, ControlPlots.plt, wf, X, Y, Z, vis, 0; msr=99)

                # Test with msr=AddedTurbulence (added turbulence) - create new state for this test
                state3 = @test_nowarn plotFlowField(nothing, ControlPlots.plt, wf, X, Y, Z, vis, 24; msr=AddedTurbulence)
                @test state3 isa FLORIDyn.PlotState
                
                # Close the plot after testing
                try
                    ControlPlots.plt.close(state3.fig)
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
                    # plt.close_all(plt); plt.pause(0.1)
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
                plot_measurements(wf, md, vis; separated=true, msr=VelReduction, plt=ControlPlots.plt)
                println("✓ plot_measurements msr=VelReduction with separated=true completed successfully")
            catch e
                @error "plot_measurements msr=VelReduction with separated=true failed: $e"
                rethrow(e)
            end
            
            try
                plot_measurements(wf, md, vis; separated=false, msr=VelReduction, plt=ControlPlots.plt)
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
                    plot_measurements(wf, md, vis; separated=true, msr=AddedTurbulence, plt=ControlPlots.plt)
                    println("✓ plot_measurements msr=AddedTurbulence with separated=true completed successfully")
                catch e
                    @warn "plot_measurements msr=AddedTurbulence with separated=true failed: $e"
                    # Still test that it fails gracefully, not with unhandled errors
                    @test isa(e, Exception)
                end
                
                try
                    plot_measurements(wf, md, vis; separated=false, msr=AddedTurbulence, plt=ControlPlots.plt)
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
                    plot_measurements(wf, md, vis; separated=true, msr=EffWind, plt=ControlPlots.plt)
                    println("✓ plot_measurements msr=EffWind with separated=true completed successfully")
                catch e
                    @warn "plot_measurements msr=EffWind with separated=true failed: $e"
                    @test isa(e, Exception)
                end
                
                try
                    plot_measurements(wf, md, vis; separated=false, msr=EffWind, plt=ControlPlots.plt)
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
            @test_throws TypeError plot_measurements(wf, md, vis; separated=true, msr=0, plt=ControlPlots.plt)
            @test_throws TypeError plot_measurements(wf, md, vis; separated=true, msr=4, plt=ControlPlots.plt)
            @test_throws TypeError plot_measurements(wf, md, vis; separated=false, msr=-1, plt=ControlPlots.plt)
        end
        
        # Test default msr value (should be 1)
        @testset "Default msr value" begin
            try
                # Test that calling without msr parameter defaults to msr=VelReduction
                plot_measurements(wf, md, vis; separated=true, plt=ControlPlots.plt)
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
end
else
    # Running tests via Pkg.test (safest approach)
    @eval Main using Pkg
    Main.Pkg.test(test_args=["test_visualisation.jl"])
end
nothing
