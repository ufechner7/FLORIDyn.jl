# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Test
using FLORIDyn
using DataFrames

# Mock plotting structure for testing
mutable struct MockPlt
    figures::Vector{String}
    methods::Dict{Symbol, Any}
    
    function MockPlt()
        new(String[], Dict{Symbol, Any}())
    end
end

function Base.getproperty(plt::MockPlt, name::Symbol)
    if name == :figures
        return getfield(plt, :figures)
    elseif name == :methods
        return getfield(plt, :methods)
    else
        # Return a generic mock function for any method
        return (args...; kwargs...) -> begin
            if name == :figure && length(args) > 0
                push!(plt.figures, string(args[1]))
            end
            return nothing
        end
    end
end

function Base.setproperty!(plt::MockPlt, name::Symbol, value)
    if name in [:figures, :methods]
        setfield!(plt, name, value)
    else
        plt.methods[name] = value
    end
end

@testset verbose=true "plotMeasurements with 54 turbines" begin
    @testset "get_layout function for large turbine counts" begin
        # Test specific cases for different turbine counts
        @test FLORIDyn.get_layout(0) == (1, 1)
        @test FLORIDyn.get_layout(1) == (1, 1)
        @test FLORIDyn.get_layout(4) == (2, 2)
        @test FLORIDyn.get_layout(6) == (2, 3)
        @test FLORIDyn.get_layout(9) == (3, 3)
        @test FLORIDyn.get_layout(12) == (3, 4)
        @test FLORIDyn.get_layout(16) == (4, 4)
        
        # Test cases for larger numbers (including 54 turbines)
        @test FLORIDyn.get_layout(20) == (4, 5)  # ceil(sqrt(20)) = 5, ceil(20/5) = 4
        @test FLORIDyn.get_layout(25) == (5, 5)  # Perfect square
        @test FLORIDyn.get_layout(54) == (7, 8)  # ceil(sqrt(54)) = 8, ceil(54/8) = 7
        @test FLORIDyn.get_layout(100) == (10, 10)  # Perfect square
        
        # Verify that the layout can accommodate all turbines
        for nT in [17, 20, 25, 30, 40, 50, 54, 60, 100]
            rows, cols = FLORIDyn.get_layout(nT)
            @test rows * cols >= nT
            @test rows > 0 && cols > 0
            # Check that it's reasonably square-like (aspect ratio not too extreme)
            aspect_ratio = max(rows, cols) / min(rows, cols)
            @test aspect_ratio <= 2.0
        end
    end
    
    @testset "plotMeasurements functionality with 54T project" begin
        # Test the core get_layout functionality with 54 turbines
        @test FLORIDyn.get_layout(54) == (7, 8)
        
        # Create a mock wind farm with 54 turbines
        wf_mock = WindFarm()
        wf_mock.nT = 54
        
        # Create mock measurement data for 54 turbines
        n_time_steps = 5
        n_turbines = 54
        
        md = DataFrame(
            Time = repeat(1.0:n_time_steps, n_turbines),
            ForeignReduction = rand(n_time_steps * n_turbines) * 0.2,
            AddedTurbulence = rand(n_time_steps * n_turbines) * 0.1,
            EffWindSpeed = rand(n_time_steps * n_turbines) * 5.0 .+ 8.0
        )
        
        # Create mock vis settings
        vis_mock = Vis()
        vis_mock.unit_test = true
        vis_mock.save = false
        vis_mock.show_plots = false
        
        # Create mock plt object with proper structure
        mock_plt = MockPlt()
        
        # Create a minimal mock ControlPlots-like object using a mutable struct
        mutable struct MockControlPlots
            plotx::Function
        end
        
        mock_controlplots = MockControlPlots((args...; kwargs...) -> begin
            push!(mock_plt.figures, "plotx_called")
            return nothing
        end)
        
        # Test VelReduction measurements with separated=true (which calls different code path)
        @test_nowarn begin
            FLORIDyn.plotMeasurements(mock_plt, wf_mock, md, vis_mock; separated=true, msr=VelReduction, pltctrl=mock_controlplots)
        end
        
        # Test AddedTurbulence measurements
        @test_nowarn begin
            FLORIDyn.plotMeasurements(mock_plt, wf_mock, md, vis_mock; separated=true, msr=AddedTurbulence, pltctrl=mock_controlplots)
        end
        
        # Test EffWind measurements  
        @test_nowarn begin
            FLORIDyn.plotMeasurements(mock_plt, wf_mock, md, vis_mock; separated=true, msr=EffWind, pltctrl=mock_controlplots)
        end
        
        # Verify that plotting functions were called
        @test length(mock_plt.figures) > 0
    end
    
    @testset "Edge cases and error handling" begin
        # Test with empty DataFrame
        empty_df = DataFrame()
        wf_mock = WindFarm()
        wf_mock.nT = 1
        vis_mock = Vis()
        vis_mock.unit_test = true
        vis_mock.save = false
        vis_mock.show_plots = false
        
        # Create mock plt
        mock_plt = Dict(
            :figure => (args...; kwargs...) -> nothing,
            :plot => (args...; kwargs...) -> nothing
        )
        
        # Test error handling for missing columns
        df_no_column = DataFrame(Time=[1.0, 2.0])
        
        @test_throws ErrorException begin
            FLORIDyn.plotMeasurements(mock_plt, wf_mock, df_no_column, vis_mock; msr=VelReduction)
        end
    end
end
