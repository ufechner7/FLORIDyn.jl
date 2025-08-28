# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test

correction = FLORIDyn.WindCorrection("None", "All", "None")
perturbation = FLORIDyn.WindPerturbation(0.0, 0.2, 0.0, 0.5, 0.0, 0.005)
shear = FLORIDyn.WindShear(0.08, 1.0)

dir_array = [
    0.0 255.0;
    20600.0 255.0;
    20900.0 195.0;
    21200.0 195.0
]

wind = FLORIDyn.Wind(
    "Constant",
    "Interpolation",
    "Constant",
    "PowerLaw",
    correction,
    perturbation,
    8.2,
    dir_array,
    0.062,
    shear
)
@testset verbose=true "correction" begin
    @testset "correctDir!" begin
        # Setup
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_All(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), false, false)
        wf = WindFarm(
            States_WF = zeros(3, 4),
            StartI = [2 2; 3 3],
            nT = 3,
            D = [126.0, 126.0, 126.0],
            nOP = 4,
            intOPs = [[1 0.5 2 0.5; 2 0.5 3 0.5]]
        )
        sim_time = 20000  # dummy input

        # Call the function
        correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)

        # Validate
        @test all(wf.States_WF[:, 2] .== 255.0)
        @test wf.States_WF[2, 4] == 255.0
    end

    @testset "correctDir! Direction_None with 4 columns" begin
        # Setup
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), false, false)
        wf = WindFarm(
            States_WF = [10.0 180.0 0.05 90.0; 8.5 190.0 0.06 100.0; 9.2 185.0 0.055 95.0],
            StartI = [2 2; 3 3],
            nT = 3,
            D = [126.0, 126.0, 126.0],
            nOP = 4,
            intOPs = [[1 0.5 2 0.5; 2 0.5 3 0.5]]
        )
        sim_time = 20000  # dummy input

        # Store initial states for comparison
        initial_direction = wf.States_WF[wf.StartI[1], 2]  # Direction at StartI position
        
        # Call the function
        correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)

        # Validate that all turbines get the same direction (phi[1])
        @test all(wf.States_WF[:, 2] .== 255.0)
        
        # Validate that OP orientation matches the turbine wind direction (not phi[1])
        # For Direction_None, OP orientation should copy the current turbine state
        @test wf.States_WF[wf.StartI[1], 4] == 255.0  # Should match updated turbine direction
    end

    @testset "correctDir! Direction_None with 3 columns" begin
        # Setup wind farm without 4th column (no OP orientation)
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), false, false)
        wf = WindFarm(
            States_WF = [10.0 180.0 0.05; 8.5 190.0 0.06; 9.2 185.0 0.055],  # Only 3 columns
            StartI = [2 2; 3 3],
            nT = 3,
            D = [126.0, 126.0, 126.0],
            nOP = 4,
            intOPs = [[1 0.5 2 0.5; 2 0.5 3 0.5]]
        )
        sim_time = 20000

        # Call the function
        correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)

        # Validate that all turbines get the same direction
        @test all(wf.States_WF[:, 2] .== 255.0)
        
        # Validate that we still have only 3 columns (no OP orientation column added)
        @test size(wf.States_WF, 2) == 3
    end

    @testset "correctDir! Direction_None vs Direction_All comparison" begin
        # Setup for Direction_None
        set_none = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), false, false)
        wf_none = WindFarm(
            States_WF = [10.0 180.0 0.05 90.0; 8.5 190.0 0.06 100.0; 9.2 185.0 0.055 95.0],
            StartI = [2 2; 3 3],
            nT = 3,
            D = [126.0, 126.0, 126.0],
            nOP = 4,
            intOPs = [[1 0.5 2 0.5; 2 0.5 3 0.5]]
        )

        # Setup for Direction_All (identical initial state)
        set_all = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_All(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), false, false)
        wf_all = WindFarm(
            States_WF = [10.0 180.0 0.05 90.0; 8.5 190.0 0.06 100.0; 9.2 185.0 0.055 95.0],
            StartI = [2 2; 3 3],
            nT = 3,
            D = [126.0, 126.0, 126.0],
            nOP = 4,
            intOPs = [[1 0.5 2 0.5; 2 0.5 3 0.5]]
        )
        sim_time = 20000

        # Apply corrections
        correctDir!(set_none.cor_dir_mode, set_none, wf_none, wind, sim_time)
        correctDir!(set_all.cor_dir_mode, set_all, wf_all, wind, sim_time)

        # Both should set all turbine directions to the same value
        @test all(wf_none.States_WF[:, 2] .== 255.0)
        @test all(wf_all.States_WF[:, 2] .== 255.0)
        @test wf_none.States_WF[:, 2] == wf_all.States_WF[:, 2]
        
        # Key difference: OP orientation handling
        # Direction_All: OP orientation = phi[1] (retrieved direction data)
        # Direction_None: OP orientation = current turbine state (after update)
        @test wf_all.States_WF[wf_all.StartI[1], 4] == 255.0  # phi[1]
        @test wf_none.States_WF[wf_none.StartI[1], 4] == 255.0  # current turbine state (same in this case)
        
        # In this test case, both should be equal since turbine state gets updated to phi[1]
        @test wf_none.States_WF[wf_none.StartI[1], 4] == wf_all.States_WF[wf_all.StartI[1], 4]
    end

    @testset "correctDir! Direction_None with different wind input modes" begin
        # Test with standard interpolation mode (already covered by other tests)
        # This test focuses on the function's robustness with different simulation times
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), false, false)
        wf = WindFarm(
            States_WF = [10.0 180.0 0.05 90.0; 8.5 190.0 0.06 100.0; 9.2 185.0 0.055 95.0],
            StartI = [2 2; 3 3],
            nT = 3,
            D = [126.0, 126.0, 126.0],
            nOP = 4,
            intOPs = [[1 0.5 2 0.5; 2 0.5 3 0.5]]
        )

        # Test with different simulation times
        for sim_time in [0, 10000, 20600, 21200]
            @test_nowarn correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)
            # Direction column should be updated
            @test all(isa.(wf.States_WF[:, 2], Float64))
            @test length(unique(wf.States_WF[:, 2])) == 1  # All turbines same direction
        end
    end

    @testset "correctDir! Direction_None return value and side effects" begin
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), false, false)
        wf = WindFarm(
            States_WF = [10.0 180.0 0.05 90.0; 8.5 190.0 0.06 100.0; 9.2 185.0 0.055 95.0],
            StartI = [2 2; 3 3],
            nT = 3,
            D = [126.0, 126.0, 126.0],
            nOP = 4,
            intOPs = [[1 0.5 2 0.5; 2 0.5 3 0.5]]
        )
        sim_time = 20000

        # Store original state for comparison (only non-direction columns)
        original_col1 = copy(wf.States_WF[:, 1])  # velocities
        original_col3 = copy(wf.States_WF[:, 3])  # turbulence

        # Function should return nothing
        result = correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)
        @test result === nothing

        # Only direction columns should be modified
        @test wf.States_WF[:, 1] == original_col1  # velocities unchanged
        @test wf.States_WF[:, 3] == original_col3  # turbulence unchanged
        @test all(wf.States_WF[:, 2] .== 255.0)    # directions updated
    end

    # @testset "correctDir! without 4th state column" begin
    #     dir_strategy = Direction_All()
    #     mock_t = MockT(zeros(3, 3), 1)
    #     wind = :mockwind
    #     sim_time = :simtime

    #     correctDir!(dir_strategy, mock_t, wind, sim_time)

    #     # @test all(mock_t.States_WF[:, 2] .== 999.9)
    #     # Should not throw an error or modify column 4
    #     # @test size(mock_t.States_WF, 2) == 3
    # end
end
nothing