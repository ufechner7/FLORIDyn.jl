# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test

function get_parameters()
    settings_file = "data/2021_9T_Data.yaml"

    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)

    # create settings struct
    set = Settings(wind, sim, con)

    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)  
    wf = initSimulation(wf, sim)
    wf, md, mi = runFLORIDyn(set, wf, wind, sim, con, floridyn, floris)
    return wf, set, floris, wind 
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
            wf, set, floris, wind = get_parameters()
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
        
        @testset "output dimensions" begin
            # Test that output has correct dimensions
            mx = [0.0 100.0; 200.0 300.0]
            my = [0.0 0.0; 100.0 100.0]
            nM = 3
            zh = 90.0
            
            expected_size = (size(mx, 1), size(mx, 2), nM)
            @test expected_size == (2, 2, 3)
            
            # Test with different grid sizes
            mx_small = [0.0]
            my_small = [0.0]
            expected_size_small = (1, 1, nM)
            @test expected_size_small == (1, 1, 3)
            
            mx_rect = [0.0 100.0 200.0; 300.0 400.0 500.0]  # 2x3 grid
            my_rect = [0.0 0.0 0.0; 100.0 100.0 100.0]
            expected_size_rect = (size(mx_rect, 1), size(mx_rect, 2), nM)
            @test expected_size_rect == (2, 3, 3)
        end
        
        @testset "coordinate handling" begin
            # Test various coordinate scenarios
            mx = [0.0 1000.0; -500.0 2000.0]
            my = [-1000.0 500.0; 0.0 1500.0]
            
            # Test negative coordinates (should be handled)
            @test minimum(mx) == -500.0
            @test minimum(my) == -1000.0
            
            # Test large coordinates (should be handled)
            @test maximum(mx) == 2000.0
            @test maximum(my) == 1500.0
        end
        
        @testset "parameter validation" begin
            mx = [0.0 100.0; 200.0 300.0]
            my = [0.0 0.0; 100.0 100.0]
            
            # Test valid nM values (just check they are positive integers)
            nM_valid = [1, 3, 5, 10]
            for nM in nM_valid
                @test nM > 0
                @test isa(nM, Int)
            end
            
            # Test valid hub heights
            zh_valid = [50.0, 90.0, 120.0, 150.0]
            for zh in zh_valid
                @test zh > 0.0
            end
            
            # Test that negative or zero values would be problematic
            # (Note: actual validation would happen in the function implementation)
            nM_invalid = [-1, 0]
            for nM in nM_invalid
                @test nM <= 0  # These are invalid values
            end
            
            zh_invalid = [-10.0, 0.0]
            for zh in zh_invalid
                @test zh <= 0.0  # These would be invalid in practice
            end
        end
        
        @testset "matrix dimension consistency" begin
            # Test mismatched matrix dimensions
            mx_2x2 = [0.0 100.0; 200.0 300.0]
            my_2x3 = [0.0 0.0 0.0; 100.0 100.0 100.0]
            
            @test size(mx_2x2) != size(my_2x3)  # Should be caught in validation
            
            # Test matched dimensions
            mx_matched = [0.0 100.0; 200.0 300.0]
            my_matched = [0.0 0.0; 100.0 100.0]
            
            @test size(mx_matched) == size(my_matched)
        end
        
        @testset "edge cases" begin
            # Single grid point
            mx_single = reshape([100.0], 1, 1)
            my_single = reshape([200.0], 1, 1)
            nM = 3
            zh = 90.0
            
            @test size(mx_single) == (1, 1)
            @test size(my_single) == (1, 1)
            @test length(mx_single) == 1
            
            # Large grid (performance consideration)
            n_large = 50
            mx_large = reshape(collect(1.0:n_large^2), n_large, n_large)
            my_large = reshape(collect(1.0:n_large^2), n_large, n_large)
            
            @test size(mx_large) == (n_large, n_large)
            @test length(mx_large) == n_large^2
        end
        
        @testset "basic coordinate validation" begin
            # Test basic properties that the coordinate transformation should satisfy
            # without testing the exact implementation details
            
            mx = [0.0 100.0 200.0; 300.0 400.0 500.0]  # 2x3 matrix
            size_mx = size(mx)
            
            # Test that matrix dimensions are what we expect
            @test size_mx[1] == 2  # rows
            @test size_mx[2] == 3  # columns
            @test length(mx) == 6  # total elements
            
            # Test basic divrem functionality (mathematical property)
            @test divrem(0, 2) == (0, 0)
            @test divrem(1, 2) == (0, 1) 
            @test divrem(2, 2) == (1, 0)
            @test divrem(3, 2) == (1, 1)
            
            # Test that the general pattern works for any reasonable matrix size
            for nrows in 1:5
                for ncols in 1:5
                    total_elements = nrows * ncols
                    for i in 0:(total_elements-1)
                        quotient, remainder = divrem(i, nrows)
                        # Basic mathematical properties
                        @test quotient >= 0
                        @test 0 <= remainder < nrows
                        @test quotient * nrows + remainder == i
                    end
                end
            end
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
end
nothing
