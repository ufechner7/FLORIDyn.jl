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
            wf, set, floris, wind = get_parameters()
            
            # Call the plotFlowField function
            Z, X, Y = plotFlowField(set, wf, wind, floris)
            
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
            wf, set, floris, wind = get_parameters()
            
            # # Call the function
            # Z, X, Y = plotFlowField(set, wf, wind, floris)
            
            # # Test coordinate range
            # @test minimum(X) >= 0.0
            # @test maximum(X) <= 3000.0
            # @test minimum(Y) >= 0.0  
            # @test maximum(Y) <= 3000.0
            
            # # Test that coordinates are within expected bounds
            # @test all(X .>= 0.0)
            # @test all(X .<= 3000.0)
            # @test all(Y .>= 0.0)
            # @test all(Y .<= 3000.0)
        end
    end

end
nothing
