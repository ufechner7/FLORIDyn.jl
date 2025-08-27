using Test
using FLORIDyn

@testset verbose=true "getDataVel function tests" begin
    
    @testset "Constant velocity mode" begin
        @testset "Basic constant velocity" begin
            settings = Settings(
                Velocity_Constant(),
                Direction_Constant(),
                TI_Constant(),
                Shear_PowerLaw(),
                Direction_Constant(),
                Velocity_Constant(),
                TI_Constant(),
                IterateOPs_basic(),
                Yaw_Constant(),
                false,
                false
            )
            
            wind_velocities = [[8.0, 9.0, 10.0]]
            turbine_indices = [1, 2]
            time_index = 1
            
            # Test function does not throw errors
            @test_nowarn getDataVel(wind_velocities, turbine_indices, time_index, settings)
            
            # Test function returns correct values for constant mode
            result = getDataVel(wind_velocities, turbine_indices, time_index, settings)
            @test result == [8.0, 9.0]  # Should return velocities for the specified turbines
        end
        
        @testset "Empty turbine array" begin
            settings = Settings(
                Velocity_Constant(),
                Direction_Constant(),
                TI_Constant(),
                Shear_PowerLaw(),
                Direction_Constant(),
                Velocity_Constant(),
                TI_Constant(),
                IterateOPs_basic(),
                Yaw_Constant(),
                false,
                false
            )
            
            wind_velocities = [[8.0, 9.0, 10.0]]
            turbine_indices = Int[]  # Empty array
            time_index = 1
            
            result = getDataVel(wind_velocities, turbine_indices, time_index, settings)
            @test result == Float64[]  # Should return empty array
        end
    end
    
    @testset "Interpolation mode" begin
        @testset "Basic interpolation" begin
            settings = Settings(
                Velocity_Interpolation(),
                Direction_Constant(),
                TI_Constant(),
                Shear_PowerLaw(),
                Direction_Constant(),
                Velocity_Constant(),
                TI_Constant(),
                IterateOPs_basic(),
                Yaw_Constant(),
                false,
                false
            )
            
            wind_velocities = [
                [8.0, 9.0, 10.0],   # Time point 1
                [8.5, 9.5, 10.5]    # Time point 2
            ]
            turbine_indices = [1, 2, 3]
            time_index = 1.5  # Non-integer for interpolation
            
            # Test function does not throw errors
            @test_nowarn getDataVel(wind_velocities, turbine_indices, time_index, settings)
            
            # Test that result is an array of correct length
            result = getDataVel(wind_velocities, turbine_indices, time_index, settings)
            @test length(result) == length(turbine_indices)
            @test isa(result, Vector{Float64})
        end
    end
    
    @testset "Return type consistency" begin
        settings = Settings(
            Velocity_Constant(),
            Direction_Constant(),
            TI_Constant(),
            Shear_PowerLaw(),
            Direction_Constant(),
            Velocity_Constant(),
            TI_Constant(),
            IterateOPs_basic(),
            Yaw_Constant(),
            false,
            false
        )
        
        # Test that function always returns Vector{Float64}
        wind_velocities = [[8.0, 9.0, 10.0]]
        
        for turbine_indices in [Int[], [1], [1, 2], [1, 2, 3]]
            result = getDataVel(wind_velocities, turbine_indices, 1, settings)
            @test isa(result, Vector{Float64})
            @test length(result) == length(turbine_indices)
        end
    end
    
    @testset "Edge cases" begin
        settings = Settings(
            Velocity_Constant(),
            Direction_Constant(),
            TI_Constant(),
            Shear_PowerLaw(),
            Direction_Constant(),
            Velocity_Constant(),
            TI_Constant(),
            IterateOPs_basic(),
            Yaw_Constant(),
            false,
            false
        )
        
        @testset "Single velocity value" begin
            wind_velocities = [[10.0]]
            turbine_indices = [1]
            time_index = 1
            
            result = getDataVel(wind_velocities, turbine_indices, time_index, settings)
            @test result == [10.0]
        end
        
        @testset "Large wind farm" begin
            wind_velocities = [collect(range(5.0, 15.0, length=100))]
            turbine_indices = collect(1:50)
            time_index = 1
            
            result = getDataVel(wind_velocities, turbine_indices, time_index, settings)
            @test length(result) == 50
            @test result == collect(range(5.0, 15.0, length=100))[1:50]
        end
    end
end
