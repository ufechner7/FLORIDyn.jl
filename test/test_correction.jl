# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test

correction = FLORIDyn.WindCorrection(vel="None", dir="All", ti="None")
perturbation = FLORIDyn.WindPerturbation(vel=false, vel_sigma=0.2, dir=false, dir_sigma=0.5, ti=false, ti_sigma=0.005)
shear = FLORIDyn.WindShear(0.08, 1.0)

dir_array = [
    0.0 255.0;
    20600.0 255.0;
    20900.0 195.0;
    21200.0 195.0
]

wind = FLORIDyn.Wind(
    input_vel="Constant",
    input_dir="Interpolation",
    input_ti="Constant",
    input_shear="PowerLaw",
    correction=correction,
    perturbation=perturbation,
    vel=8.2,
    dir=dir_array,
    ti=0.062,
    shear=shear
)
sim_time = 20600.0  # default simulation time for initial direction tests
@testset verbose=true "correction" begin
    @testset "correctDir!" begin
        # Setup
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_All(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
        wf = WindFarm(
            nT = 3,
            nOP = 3,
            States_WF = [8.0 200.0 0.05 0.0; 9.0 205.0 0.06 0.0; 7.5 210.0 0.07 0.0],
            StartI = [1 1; 2 2; 3 3],
            D = [120.0,120.0,120.0],
            intOPs = [ Matrix{Float64}(undef,0,0), [1.0 0.4 2.0 0.6; 2.0 0.3 3.0 0.7], Matrix{Float64}(undef,0,0) ],
            dep = [ Int[], [1,2], Int[] ],
            Weight = [ Float64[], [0.2, 0.8], Float64[] ],
        )
        # Call the function
        correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)

        # Validate
        @test all(wf.States_WF[:, 2] .== 255.0)
        @test wf.States_WF[2, 4] == 255.0
    end

    @testset "correctDir! Direction_None with 4 columns" begin
        # Setup
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                        Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
        wf = WindFarm(
            nT = 2,
            nOP = 2,
            States_WF = [8.0 200.0 0.05 0.0; 9.0 205.0 0.06 0.0],
            StartI = [1 1; 2 2],
            D = [120.0,120.0],
            intOPs = [ Matrix{Float64}(undef,0,0), [1.0 0.5 2.0 0.5] ],
            dep = [ Int[], [1,2] ],
            Weight = [ Float64[], [0.0] ],
        )
        # Call the function
        correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)

        # Validate that only starting positions get the same direction (phi[1])
        # Direction_None only updates starting positions, not all operational points
        @test wf.States_WF[wf.StartI[1, 1], 2] == 255.0  # First turbine starting position
        @test wf.States_WF[wf.StartI[2, 1], 2] == 255.0  # Second turbine starting position
        # For Direction_None, OP orientation should be set to phi[1] directly
        @test wf.States_WF[wf.StartI[1, 1], 4] == 255.0
        @test wf.States_WF[wf.StartI[2, 1], 4] == 255.0
    end

    @testset "correctDir! Direction_None with 3 columns" begin
        # Setup wind farm without 4th column (no OP orientation)
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
  
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

        # Validate that only starting positions get the direction update
        # Direction_None only updates starting positions [2, 3], not all rows
        @test wf.States_WF[2, 2] == 255.0  # Starting position for first turbine
        @test wf.States_WF[3, 2] == 255.0  # Starting position for second turbine
        @test wf.States_WF[1, 2] == 180.0  # Non-starting position should remain unchanged
        
        # Validate that we still have only 3 columns (no OP orientation column added)
        @test size(wf.States_WF, 2) == 3
    end

    @testset "correctDir! Direction_None vs Direction_All comparison" begin
        # Setup for Direction_None
        set_none = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
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
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
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

        # Key difference: Direction_None only updates starting positions, Direction_All updates all OPs
        # Direction_None: Only starting positions [2,3] should be updated to 255.0
        @test wf_none.States_WF[2, 2] == 255.0  # Starting position updated
        @test wf_none.States_WF[3, 2] == 255.0  # Starting position updated  
        @test wf_none.States_WF[1, 2] == 180.0  # Non-starting position unchanged
        
        # Direction_All: All operational points should be updated to 255.0
        @test all(wf_all.States_WF[:, 2] .== 255.0)
        
        # The arrays should NOT be equal due to different update strategies
        @test wf_none.States_WF[:, 2] != wf_all.States_WF[:, 2]
        
        # OP orientation handling
        # Direction_All: OP orientation = updated turbine state for starting positions
        # Direction_None: OP orientation = phi[1] directly for starting positions
        @test wf_all.States_WF[wf_all.StartI[1, 1], 4] == wf_all.States_WF[wf_all.StartI[1, 1], 2]  # Copies turbine state
        @test wf_none.States_WF[wf_none.StartI[1, 1], 4] == 255.0  # Set to phi[1] directly
        
        # In this case, both OP orientations should be equal since both end up as 255.0
        @test wf_none.States_WF[wf_none.StartI[1, 1], 4] == wf_all.States_WF[wf_all.StartI[1, 1], 4]
    end

    @testset "correctDir! Direction_None with different wind input modes" begin
        # Test with standard interpolation mode (already covered by other tests)
        # This test focuses on the function's robustness with different simulation times
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
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
            # Direction_None only updates starting positions, so not all turbines have same direction
            # Check that starting positions are updated correctly
            @test all(isa.(wf.States_WF[wf.StartI, 2], Float64))
            # Starting positions should have same direction value
            starting_directions = wf.States_WF[wf.StartI, 2]
            @test length(unique(starting_directions)) == 1  # Starting positions same direction
        end
    end

    @testset "correctDir! Direction_None return value and side effects" begin
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_None(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
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
        original_col2 = copy(wf.States_WF[:, 2])  # directions (to check selective update)

        # Function should return nothing
        result = correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)
        @test result === nothing

        # Only direction columns should be modified, and only for starting positions
        @test wf.States_WF[:, 1] == original_col1  # velocities unchanged
        @test wf.States_WF[:, 3] == original_col3  # turbulence unchanged
        
        # Direction_None only updates starting positions
        @test wf.States_WF[2, 2] == 255.0  # Starting position updated
        @test wf.States_WF[3, 2] == 255.0  # Starting position updated
        @test wf.States_WF[1, 2] == original_col2[1]  # Non-starting position unchanged
    end

    @testset "correctDir! Direction_Influence basic cases" begin
        # Reuse global wind (dir_array) so phi at chosen times is 255.0
        set_infl = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_Influence(), 
                        Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
        sim_time = 50.0  # Within first plateau => phi = 255

        # Case 1: No dependencies -> raw phi applied
        wf1 = WindFarm(
            nT = 2,
            nOP = 4,
            States_WF = [0.0 180.0 0.05 0.0; 0.0 170.0 0.05 0.0],
            StartI = [1 2],  # 1 x nT layout (each turbine's start row)
            D = [120.0, 120.0],
            intOPs = Matrix{Float64}[],
            Weight = Vector{Float64}[],
            dep = [Int[], Int[]],
        )
        correctDir!(set_infl.cor_dir_mode, set_infl, wf1, wind, sim_time)
        @test all(wf1.States_WF[:, 2] .== 255.0)
        @test all(wf1.States_WF[:, 4] .== 255.0)

        # Base pre-update states for later cases
        base_states = [0.0 260.0 0.05 0.0; 0.0 200.0 0.05 0.0; 0.0 210.0 0.05 0.0]

        # Case 2: Single intOP row (two indices + weights) for turbine 2
        wf2 = WindFarm(
            nT = 3,
            nOP = 6,
            States_WF = copy(base_states),
            StartI = [1 2 3],
            D = [120.0, 120.0, 120.0],
            intOPs = [Matrix{Float64}(undef,0,0), [1 0.25 3 0.75], Matrix{Float64}(undef,0,0)],
            Weight = [Float64[], Float64[], Float64[]],
            dep = [Int[], [1,3], Int[]],
        )
        correctDir!(set_infl.cor_dir_mode, set_infl, wf2, wind, sim_time)
        expected_dir = 255.0*0.25 + 210.0*0.75  # turbine 1 updated to 255 before combining
        @test isapprox(wf2.States_WF[2, 2], expected_dir; atol=1e-8)
        @test wf2.States_WF[1,2] == 255.0
        @test wf2.States_WF[3,2] == 255.0

        # Case 3: Multiple intOP rows with non-zero weights
        int_rows = [1 0.5 2 0.5; 1 0.8 2 0.2]
        weights = [0.6, 0.4]
        wf3 = WindFarm(
            nT = 3,
            nOP = 6,
            States_WF = [0.0 260.0 0.05 0.0; 0.0 200.0 0.05 0.0; 0.0 180.0 0.05 0.0],
            StartI = [1 2 3],
            D = [120.0, 120.0, 120.0],
            intOPs = [Matrix{Float64}(undef,0,0), Matrix{Float64}(undef,0,0), int_rows],
            Weight = [Float64[], Float64[], weights],
            dep = [Int[], Int[], [1,2]],
        )
        correctDir!(set_infl.cor_dir_mode, set_infl, wf3, wind, sim_time)
        @test wf3.States_WF[1,2] == 255.0
        @test wf3.States_WF[2,2] == 255.0
        local1 = 255.0*0.5 + 255.0*0.5
        local2 = 255.0*0.8 + 255.0*0.2
        expected3 = (0.6*local1 + 0.4*local2)/(0.6+0.4)
        @test isapprox(wf3.States_WF[3,2], expected3; atol=1e-8)

        # Case 4: Multiple rows but zero weights -> fallback to phi
        wf4 = WindFarm(
            nT = 2,
            nOP = 4,
            States_WF = [0.0 260.0 0.05 0.0; 0.0 200.0 0.05 0.0],
            StartI = [1 2],
            D = [120.0, 120.0],
            intOPs = [Matrix{Float64}(undef,0,0), [1 0.5 2 0.5; 1 0.5 2 0.5]],
            Weight = [Float64[], [0.0, 0.0]],
            dep = [Int[], [1,2]],
        )
        correctDir!(set_infl.cor_dir_mode, set_infl, wf4, wind, sim_time)
        @test wf4.States_WF[2,2] == 255.0

        # Case 5: Malformed intOPs shape -> fallback to phi
        wf5 = WindFarm(
            nT = 2,
            nOP = 4,
            States_WF = [0.0 260.0 0.05 0.0; 0.0 200.0 0.05 0.0],
            StartI = [1 2],
            D = [120.0, 120.0],
            intOPs = [Matrix{Float64}(undef,0,0), [1 0.5 2]],  # wrong shape
            Weight = [Float64[], Float64[]],
            dep = [Int[], [1,2]],
        )
        correctDir!(set_infl.cor_dir_mode, set_infl, wf5, wind, sim_time)
        @test wf5.States_WF[2,2] == 255.0

        # Case 6: Orientation only updated at StartI rows (additional non-start rows present)
        wf6 = WindFarm(
            nT = 2,
            nOP = 4,
            States_WF = [0.0 180.0 0.05 1.0; 0.0 190.0 0.05 2.0; 0.0 200.0 0.05 3.0; 0.0 210.0 0.05 4.0],
            StartI = [1 3],  # turbine 1 -> row 1, turbine 2 -> row 3
            D = [120.0, 120.0],
            intOPs = Matrix{Float64}[],
            Weight = Vector{Float64}[],
            dep = [Int[], Int[]],
        )
        correctDir!(set_infl.cor_dir_mode, set_infl, wf6, wind, sim_time)
        # Directions only at start rows updated to phi
        @test wf6.States_WF[1,2] == 255.0
        @test wf6.States_WF[3,2] == 255.0
        @test wf6.States_WF[2,2] == 190.0
        @test wf6.States_WF[4,2] == 210.0
        # Orientation updated only at StartI rows
        @test wf6.States_WF[1,4] == 255.0
        @test wf6.States_WF[3,4] == 255.0
        @test wf6.States_WF[2,4] == 2.0
        @test wf6.States_WF[4,4] == 4.0
    end

    @testset "correctDir! Direction_All variants" begin
        # 4-column case
        set_all = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_All(), 
                        Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
        wf_a = WindFarm(
            States_WF = [10.0 180.0 0.05 90.0; 8.5 200.0 0.06 100.0; 9.2 220.0 0.055 95.0],
            StartI = [1 1; 2 2; 3 3],
            nT = 3,
            D = [126.0,126.0,126.0],
            nOP = 4,
            intOPs = Matrix{Float64}[]
        )
        correctDir!(set_all.cor_dir_mode, set_all, wf_a, wind, 20600.0)
        @test all(wf_a.States_WF[:,2] .== 255.0)
        @test all(wf_a.States_WF[wf_a.StartI,4] .== 255.0)

        # 3-column case (no orientation column)
        wf_b = WindFarm(
            States_WF = [10.0 180.0 0.05; 8.5 200.0 0.06; 9.2 220.0 0.055],
            StartI = [1 1; 2 2; 3 3],
            nT = 3,
            D = [126.0,126.0,126.0],
            nOP = 4,
            intOPs = Matrix{Float64}[]
        )
        correctDir!(set_all.cor_dir_mode, set_all, wf_b, wind, 20600.0)
        @test all(wf_b.States_WF[:,2] .== 255.0)
        @test size(wf_b.States_WF,2) == 3  # unchanged column count

        # Interpolation over time (verify varying phi applied uniformly)
        wf_c = WindFarm(
            States_WF = [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0],
            StartI = [1 1; 2 2],
            nT = 2,
            D = [120.0,120.0],
            nOP = 2,
            intOPs = Matrix{Float64}[]
        )
        # Before change segment end
        correctDir!(set_all.cor_dir_mode, set_all, wf_c, wind, 20600.0)
        @test all(wf_c.States_WF[:,2] .== 255.0)
        # Mid-transition between 20600 and 20900 should linearly interpolate (halfway at 20750)
        correctDir!(set_all.cor_dir_mode, set_all, wf_c, wind, 20750.0)
        @test all(isapprox.(wf_c.States_WF[:,2], 225.0; atol=1e-6))
        # After transition end
        correctDir!(set_all.cor_dir_mode, set_all, wf_c, wind, 20900.0)
        @test all(wf_c.States_WF[:,2] .== 195.0)

        # Side effects: only direction (and orientation if present) changes
        wf_d = WindFarm(
            States_WF = [10.0 180.0 0.05 90.0; 8.5 200.0 0.06 100.0],
            StartI = [1 1; 2 2],
            nT = 2,
            D = [126.0,126.0],
            nOP = 2,
            intOPs = Matrix{Float64}[]
        )
        vel_before = copy(wf_d.States_WF[:,1])
        turb_before = copy(wf_d.States_WF[:,3])
        res = correctDir!(set_all.cor_dir_mode, set_all, wf_d, wind, 0.0)
        @test res === nothing
        @test wf_d.States_WF[:,1] == vel_before
        @test wf_d.States_WF[:,3] == turb_before
        @test all(wf_d.States_WF[:,2] .== 255.0)
        @test all(wf_d.States_WF[wf_d.StartI,4] .== 255.0)
    end

    @testset "correctTI! TI_None basic cases" begin
        # Use TI interpolation so we exercise getWindTiT path
        TI_array = [0.0 0.075; 100.0 0.075]
        wind_ti_none = FLORIDyn.Wind(
            input_vel = "Constant",
            input_dir = "Interpolation",
            input_ti = "Interpolation",
            input_shear = "PowerLaw",
            correction = correction,
            perturbation = perturbation,
            vel = 8.2,
            dir = [0.0 250.0; 10.0 250.0],
            ti = TI_array,
            shear = shear,
        )
        set_ti_none = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Interpolation(), Shear_PowerLaw(), Direction_None(), Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
        t = 50.0

        # Case 1: 4-column States_WF, TI column initially zero
        wf_a_ti = WindFarm(
            nT = 2,
            nOP = 2,
            States_WF = [10.0 180.0 0.0 90.0; 9.5 185.0 0.0 95.0],
            StartI = [1 1; 2 2],
            D = [120.0,120.0],
            intOPs = Matrix{Float64}[],
            dep = [Int[], Int[]],
        )
        correctTI!(set_ti_none.cor_turb_mode, set_ti_none, wf_a_ti, wind_ti_none, t)
        @test all(isapprox.(wf_a_ti.States_WF[:,3], 0.075; atol=1e-8))
        # Other columns unchanged
        @test wf_a_ti.States_WF[:,1] == [10.0, 9.5]
        @test wf_a_ti.States_WF[:,2] == [180.0, 185.0]

        # Case 2: 3-column States_WF (no orientation column)
        wf_b_ti = WindFarm(
            nT = 3,
            nOP = 3,
            States_WF = [8.0 170.0 0.01; 8.2 172.0 0.02; 8.1 171.0 0.03],
            StartI = [1 1; 2 2; 3 3],
            D = [120.0,120.0,120.0],
            intOPs = Matrix{Float64}[],
            dep = [Int[], Int[], Int[]],
        )
        correctTI!(set_ti_none.cor_turb_mode, set_ti_none, wf_b_ti, wind_ti_none, t)
        @test all(isapprox.(wf_b_ti.States_WF[:,3], 0.075; atol=1e-8))
        @test size(wf_b_ti.States_WF,2) == 3

        # Case 3: Return value is nothing
        ret = correctTI!(set_ti_none.cor_turb_mode, set_ti_none, wf_b_ti, wind_ti_none, t)
        @test ret === nothing
    end

    @testset "correctTI! TI_Influence basic cases" begin
        # Turbulence interpolation data (constant 0.10 across time)
        TI_array = [0.0 0.10; 100.0 0.10]
        wind_ti = FLORIDyn.Wind(
            input_vel = "Constant",
            input_dir = "Interpolation",
            input_ti = "Interpolation",
            input_shear = "PowerLaw",
            correction = correction,
            perturbation = perturbation,
            vel = 8.2,
            dir = [0.0 250.0; 10.0 250.0],
            ti = TI_array,
            shear = shear,
        )
        set_infl_ti = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Interpolation(), Shear_PowerLaw(), Direction_None(), Velocity_None(), TI_Influence(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
        t = 50.0

        # Case 1: No dependencies
        wf1 = WindFarm(
            nT = 2,
            nOP = 2,
            States_WF = [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0],
            StartI = [1 2],  # 1×nT layout
            D = [120.0, 120.0],
            intOPs = Matrix{Float64}[],
            dep = [Int[], Int[]],
        )
        correctTI!(set_infl_ti.cor_turb_mode, set_infl_ti, wf1, wind_ti, t)
        @test all(isapprox.(wf1.States_WF[:,3], 0.10; atol=1e-8))

        # Seed upstream TI values for influence tests
        wf_base_states = [0.0 0.0 0.08 0.0; 0.0 0.0 0.12 0.0; 0.0 0.0 0.16 0.0]

        # Case 2: Single intOP row
        wf2 = WindFarm(
            nT = 3,
            nOP = 6,
            States_WF = copy(wf_base_states),
            StartI = [1 2 3],  # 1×nT layout
            D = [120.0,120.0,120.0],
            intOPs = [Matrix{Float64}(undef,0,0), [1 0.25 3 0.75], Matrix{Float64}(undef,0,0)],
            dep = [Int[], [1,3], Int[]],
        )
        correctTI!(set_infl_ti.cor_turb_mode, set_infl_ti, wf2, wind_ti, t)
        expected_ti = 0.10*0.25 + 0.16*0.75  # turbine 1 updated to 0.10 first
        @test isapprox(wf2.States_WF[2,3], expected_ti; atol=1e-8)

        # Case 3: Multiple rows -> mean of per-row interpolations
        int_rows = [1 0.5 2 0.5; 1 0.2 2 0.8]
        wf3 = WindFarm(
            nT = 3,
            nOP = 6,
            States_WF = copy(wf_base_states),
            StartI = [1 2 3],  # 1×nT layout
            D = [120.0,120.0,120.0],
            intOPs = [Matrix{Float64}(undef,0,0), Matrix{Float64}(undef,0,0), int_rows],
            dep = [Int[], Int[], [1,2]],
        )
        correctTI!(set_infl_ti.cor_turb_mode, set_infl_ti, wf3, wind_ti, t)
        local1 = 0.10*0.5 + 0.10*0.5
        local2 = 0.10*0.2 + 0.10*0.8
        expected3 = (local1 + local2)/2
        @test isapprox(wf3.States_WF[3,3], expected3; atol=1e-8)

        # Case 4: Malformed intOPs -> fallback to base TI
        wf4 = WindFarm(
            nT = 2,
            nOP = 4,
            States_WF = [0.0 0.0 0.07 0.0; 0.0 0.0 0.09 0.0],
            StartI = [1 2],  # 1×nT layout
            D = [120.0,120.0],
            intOPs = [Matrix{Float64}(undef,0,0), [1 0.5 2]],
            dep = [Int[], [1,2]],
        )
        correctTI!(set_infl_ti.cor_turb_mode, set_infl_ti, wf4, wind_ti, t)
        @test isapprox(wf4.States_WF[2,3], 0.10; atol=1e-8)
    end

    @testset "correctVel Velocity_Influence basic cases" begin
        correction_vi = FLORIDyn.WindCorrection("Constant","Constant","Constant")
        perturbation_vi = FLORIDyn.WindPerturbation(false,0.0,false,0.0,false,0.0)
        set_vi = Settings(Velocity_Constant(), Direction_Constant(), TI_Constant(), Shear_PowerLaw(), Direction_None(), Velocity_Influence(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
        floris_dummy = FLORIDyn.Floris(
            alpha = 1.88,
            beta = 0.0,
            k_a = 0.0,
            k_b = 0.0,
            k_fa = 0.0,
            k_fb = 0.0,
            k_fc = 0.0,
            k_fd = 0.0,
            eta = 1,
            p_p = 1.0,
            airDen = 1.225,
            TIexp = 1,
            rotor_points = nothing,
        )
        tmpM = ones(3,1)
        wind_vi = FLORIDyn.Wind(
            input_vel = "Constant",
            input_dir = "Constant",
            input_ti = "Constant",
            input_shear = "PowerLaw",
            correction = correction_vi,
            perturbation = perturbation_vi,
            vel = 8.0,
            dir = [0.0 250.0],
            ti = 0.08,
            shear = nothing,
        )

        # Case 1: No dependencies
        wf_v1 = WindFarm(
            nT = 2,
            nOP = 2,
            States_WF = [5.0 200.0 0.05; 6.0 205.0 0.06],
            StartI = [1 2],  # 1×nT layout
            D = [120.0,120.0],
            intOPs = Matrix{Float64}[],
            dep = [Int[], Int[]],
            Weight = Vector{Float64}[],
        )
    FLORIDyn.correctVel!(Velocity_Influence(), set_vi, wf_v1, wind_vi, 10.0, floris_dummy, tmpM)
        @test all(isapprox.(wf_v1.States_WF[:,1], 8.0; atol=1e-8))

        # Case 2: Single-row interpolation for turbine 1
        wf_v2 = WindFarm(
            nT = 2,
            nOP = 2,
            States_WF = [8.0 200.0 0.05; 7.0 205.0 0.06],
            StartI = [1 2],  # 1×nT layout
            D = [120.0,120.0],
            intOPs = [ [1.0 0.25 2.0 0.75], Matrix{Float64}(undef,0,0) ],
            dep = [ [1,2], Int[] ],
            Weight = Vector{Float64}[],
        )
        # Expected uses pre-correction states
        pre2 = copy(wf_v2.States_WF[:,1])
        exp2 = pre2[1]*0.25 + pre2[2]*0.75
        FLORIDyn.correctVel!(Velocity_Influence(), set_vi, wf_v2, wind_vi, 10.0, floris_dummy, tmpM)
        @test isapprox(wf_v2.States_WF[1,1], exp2; atol=1e-8)

        # Case 3: Multi-row with weights
        wf_v3 = WindFarm(
            nT = 3,
            nOP = 3,
            States_WF = [8.0 200.0 0.05; 9.0 205.0 0.06; 7.5 210.0 0.07],
            StartI = [1 2 3],  # 1×nT layout
            D = [120.0,120.0,120.0],
            intOPs = [ Matrix{Float64}(undef,0,0), [1.0 0.4 2.0 0.6; 2.0 0.3 3.0 0.7], Matrix{Float64}(undef,0,0) ],
            dep = [ Int[], [1,2], Int[] ],
            Weight = [ Float64[], [0.2, 0.8], Float64[] ],
        )
    # Compute expected before correction (turbine 1 has no deps -> will become base 8.0 first)
    pre3 = copy(wf_v3.States_WF[:,1])
    row1 = 8.0*0.4 + pre3[2]*0.6            # indices 1 (after correction) and 2 (pre-correction)
    row2 = pre3[2]*0.3 + pre3[3]*0.7        # indices 2 and 3 pre-correction
    exp3 = (0.2*row1 + 0.8*row2)/(0.2+0.8)
    FLORIDyn.correctVel!(Velocity_Influence(), set_vi, wf_v3, wind_vi, 10.0, floris_dummy, tmpM)
    @test isapprox(wf_v3.States_WF[2,1], exp3; atol=1e-8)

        # Case 4: Zero weights fallback
        wf_v4 = WindFarm(
            nT = 2,
            nOP = 2,
            States_WF = [8.0 200.0 0.05; 9.0 205.0 0.06],
            StartI = [1 2],  # 1×nT layout
            D = [120.0,120.0],
            intOPs = [ Matrix{Float64}(undef,0,0), [1.0 0.5 2.0 0.5] ],
            dep = [ Int[], [1,2] ],
            Weight = [ Float64[], [0.0] ],
        )
        FLORIDyn.correctVel!(Velocity_Influence(), set_vi, wf_v4, wind_vi, 10.0, floris_dummy, tmpM)
        @test isapprox(wf_v4.States_WF[2,1], 8.0; atol=1e-8)
    end

    @testset "correctVel Velocity_None basic cases" begin
        # Settings and inputs (constant velocity 9.5 m/s)
        set_vn = Settings(Velocity_Constant(), Direction_Constant(), TI_Constant(), Shear_PowerLaw(), Direction_None(), Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA(), Induction_Constant(), false, false)
        correction_vn = FLORIDyn.WindCorrection("Constant","Constant","Constant")
        perturbation_vn = FLORIDyn.WindPerturbation(false,0.0,false,0.0,false,0.0)
        wind_vn = FLORIDyn.Wind(
            input_vel = "Constant",
            input_dir = "Constant",
            input_ti = "Constant",
            input_shear = "PowerLaw",
            correction = correction_vn,
            perturbation = perturbation_vn,
            vel = 9.5,
            dir = [0.0 270.0],
            ti = 0.09,
            shear = nothing,
        )
        floris_dummy = FLORIDyn.Floris(
            alpha = 1.88, beta = 0.0,
            k_a = 0.0, k_b = 0.0, k_fa = 0.0, k_fb = 0.0, k_fc = 0.0, k_fd = 0.0,
            eta = 1, p_p = 1.0, airDen = 1.225, TIexp = 1, rotor_points = nothing,
        )
        tmpM = ones(2,1)

        # Case 1: 4-column States_WF (velocity, direction, ti, orientation)
        wf4 = WindFarm(
            nT = 2,
            nOP = 2,
            States_WF = [5.0 200.0 0.05 10.0; 6.0 205.0 0.06 20.0],
            StartI = [1 2],   # treat as vector of indices so wf.States_WF[StartI,1] is length 2
            D = [120.0,120.0],
            intOPs = Matrix{Float64}[],
        )
        ret_wf4, ret_wind4 = FLORIDyn.correctVel!(Velocity_None(), set_vn, wf4, wind_vn, 100.0, floris_dummy, tmpM)
        @test ret_wf4 === wf4
        @test ret_wind4 === wind_vn
        @test all(isapprox.(wf4.States_WF[:,1], 9.5; atol=1e-8))
        # Non-velocity columns unchanged
        @test wf4.States_WF[:,2] == [200.0,205.0]
        @test wf4.States_WF[:,3] == [0.05,0.06]
        @test wf4.States_WF[:,4] == [10.0,20.0]

        # Case 2: 3-column States_WF (no orientation column)
        wf3 = WindFarm(
            nT = 3,
            nOP = 3,
            States_WF = [7.0 180.0 0.07; 8.0 190.0 0.08; 6.5 185.0 0.09],
            StartI = [1 2 3],  # vector form
            D = [120.0,120.0,120.0],
            intOPs = Matrix{Float64}[],
        )
        FLORIDyn.correctVel!(Velocity_None(), set_vn, wf3, wind_vn, 250.0, floris_dummy, ones(3,1))
        @test all(isapprox.(wf3.States_WF[:,1], 9.5; atol=1e-8))
        @test size(wf3.States_WF,2) == 3  # still 3 columns
        @test wf3.States_WF[:,2] == [180.0,190.0,185.0]
        @test wf3.States_WF[:,3] == [0.07,0.08,0.09]
    end
end
nothing