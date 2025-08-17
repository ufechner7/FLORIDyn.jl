# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause
 
using FLORIDyn, Test, Statistics, MAT, DataFrames

matlab_file   = "test/data/input_setUpTmpWFAndRun_2_steps.mat"
after_interpolateOPs_T_file = "test/data/after_interpolateOPs_T.mat"
vars = matread(matlab_file)
vars_after_interpolateOPs_T = matread(after_interpolateOPs_T_file)
wf_dict = vars["T"]
wf_dict_ref = vars_after_interpolateOPs_T["T"]

if !isdefined(Main, :TestHelpers)
    include("test_helpers.jl")
end
using .TestHelpers

@testset verbose = true "Setup and Run Tests" begin
    @testset "compare_windFarms_function" begin
        # Test with the converted WindFarm and a copy
        wf_converted = wf_dict2windfarm(wf_dict)
        wf_copy = deepcopy(wf_converted)
        
        # Test identical WindFarms
        @test compare_windFarms(wf_converted, wf_copy, detailed=false) == true
        
        # Test with slight modification
        wf_modified = deepcopy(wf_converted)
        wf_modified.nT = 10  # Change number of turbines
        
        @test compare_windFarms(wf_converted, wf_modified, detailed=false) == false
        
        # Test with floating-point difference within tolerance
        wf_float_diff = deepcopy(wf_converted)
        if !isempty(wf_float_diff.States_WF)
            wf_float_diff.States_WF[1,1] += 1e-12  # Very small difference
        end
        
        @test compare_windFarms(wf_converted, wf_float_diff, detailed=false, tolerance=1e-10) == true
        @test compare_windFarms(wf_converted, wf_float_diff, detailed=false, tolerance=1e-14) == false
        
        # Test with name differences
        wf_name_diff = deepcopy(wf_converted)
        if !isempty(wf_name_diff.Names_T)
            wf_name_diff.Names_T[1] = "different_name"
        end
        
        @test compare_windFarms(wf_converted, wf_name_diff, detailed=false) == false
        
        # Test detailed output (just verify it runs without error)
        # println("\nTesting detailed comparison output:")
        # compare_windFarms(wf_converted, wf_copy, detailed=true)
    end

    @testset "wf_dict_to_windfarm_conversion" begin
        # Test the conversion function with MAT file data
        wf_converted = wf_dict2windfarm(wf_dict)
        
        @test wf_converted isa WindFarm
        @test wf_converted.nT == 9
        @test wf_converted.nOP == 200
        @test wf_converted.Names_T == ["a", "yaw", "TI"]
        @test wf_converted.Names_WF == ["wind_vel", "wind_dir", "TI0", "OP_ori"] 
        @test wf_converted.Names_OP == ["x0", "y0", "z0", "x1", "y1", "z1"]
        @test size(wf_converted.posBase) == (9, 3)
        @test length(wf_converted.D) == 9
        @test length(wf_converted.intOPs) == 9
        @test length(wf_converted.Weight) == 9
        @test length(wf_converted.dep) == 9
        
        # Test DataFrame property access
        turbines_df = wf_converted.turbines
        @test turbines_df isa DataFrame
        @test size(turbines_df) == (1800, 5)  # nOP×nT = 200×9 = 1800 rows, 5 columns (OP + Turbine + 3 states)
        @test names(turbines_df) == ["OP", "Turbine", "a", "yaw", "TI"]
        
        wf_df = wf_converted.windfield
        @test wf_df isa DataFrame
        @test size(wf_df) == (1800, 4)  # 1800 timesteps, 4 wind field variables
        
        ops_df = wf_converted.ops
        @test ops_df isa DataFrame
        @test size(ops_df) == (1800, 6)  # 1800 states, 6 OP variables
    end
    
    @testset "setUpTmpWFAndRun_basic" begin
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn, ta = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        wf.dep = findTurbineGroups(wf, floridyn)
        # Create unified buffers for interpolateOPs!
        intOPs_buffers = [zeros(length(wf.dep[iT]), 4) for iT in 1:wf.nT]
        unified_buffers = create_unified_buffers(wf)
        wf.intOPs = interpolateOPs!(unified_buffers, intOPs_buffers, wf)
        wf_old = deepcopy(wf)
        M, wf = setUpTmpWFAndRun(set, wf, floris, wind)
        @test ! structs_equal(wf_old, wf; prn=false)
    end

    @testset "setUpTmpWFAndRun_vs_matlab" begin
        settings_file = "data/2021_9T_Data.yaml"
        # get the settings for the wind field, simulator and controller
        wind, sim, con, floris, floridyn, ta = setup(settings_file)
        # create settings struct
        set = Settings(wind, sim, con)
        wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
        sim.n_sim_steps = 2
        wf1 = wf_dict2windfarm(wf_dict)
        tmpM, wf = setUpTmpWFAndRun(set, wf1, floris, wind)
        tmpM_ref = [0, 9.094956525547340e-02, 0.1184213787906785, 0, 9.094956525547351e-02,
                    0.1184213787906786, 0, 9.094956525547340e-02, 1.010889141095400e-01]
        @test all(tmpM[:,2] .≈ tmpM_ref)
    end
end
nothing
