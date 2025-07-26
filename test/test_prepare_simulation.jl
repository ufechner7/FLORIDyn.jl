# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test, LinearAlgebra



@testset "prepare_simulation                                      " begin
    settings_file = "data/2021_9T_Data.yaml"
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn = setup(settings_file)
    # create settings struct
    set = Settings(wind, sim, con)
    # % Load linked data
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)
    @test wind.vel == 8.2
    @test wind.dir ≈ [
                        0.0    255.0
                    20600.0  255.0
                    20900.0  195.0
                    21200.0  195.0
                    ]
    @test wind.ti == 0.0620
    @test wind.shear.z0 == 1.0
    @test wind.shear.alpha == 0.08
    # Define the expected matrix
    expected_posBase = [
        600.0  2400.0    0.0;
       1500.0  2400.0    0.0;
       2400.0  2400.0    0.0;
        600.0  1500.0    0.0;
       1500.0  1500.0    0.0;
       2400.0  1500.0    0.0;
        600.0   600.0    0.0;
       1500.0   600.0    0.0;
       2400.0   600.0    0.0
    ]

    @test wf.posBase == expected_posBase
    @test wf.nT == 9
    expected = [
        0.0 0.0 119.0
        0.0 0.0 119.0
        0.0 0.0 119.0
        0.0 0.0 119.0
        0.0 0.0 119.0
        0.0 0.0 119.0
        0.0 0.0 119.0
        0.0 0.0 119.0
        0.0 0.0 119.0
    ]

    @test wf.posNac    == expected

    @test wf.D         == fill(178.4, 9)
    @test wf.Names_OP  == ["x0", "y0", "z0", "x1", "y1", "z1"]
    @test wf.Names_T   == ["a", "yaw", "TI"]
    @test wf.Names_WF  == ["wind_vel", "wind_dir", "TI0", "OP_ori"]
    @test wf.StartI    == [1 201 401 601 801 1001 1201 1401 1601]
    @test wf.nOP       == 200
    @test wf.red_arr   == ones(9,9)

    @test floridyn.deltaUW == 10.0
    @test size(con.yaw_data) == (603, 10)
    @test sum(con.yaw_data) ≈ 1.37333255e7
    @test con.tanh_yaw == false

    @test size(wf.States_OP) == (1800, 6)
    @test sum(wf.States_OP) ≈ 1.8683e+07 rtol=1e-4

    @test size(wf.States_T) == (1800, 3)
    @test sum(wf.States_T) ≈ 702 rtol=1e-4

    @test size(wf.States_WF) == (1800, 4)
    @test sum(wf.States_WF) ≈ 9.3287e+05 rtol=1e-4

    @test sim.n_sim_steps == 301
    @test floris.rotor_points == 50

    # Test with constant wind direction
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_dir = "Constant"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    # Test with InterpTurbine
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_dir = "InterpTurbine"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    # Test with constant Interpolation_wErrorCov
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_dir = "Interpolation_wErrorCov"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)
    
    # Test with Constant_wErrorCov
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_dir = "Constant_wErrorCov"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    # Test with input_ti Interpolation
    # This will generate a demo CSV file if it does not exist
    rm("data/2021_9T_Data/WindTI.csv", force=true)
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_ti = "Interpolation"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    #  wind.input_dir == "RW_with_Mean"
    # TODO: Compare the result with the Matlab version
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_dir = "RW_with_Mean"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    # wind.input_dir == "InterpTurbine_wErrorCov"
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_dir = "InterpTurbine_wErrorCov"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    # input_ti == "InterpTurbine"
    # TODO: Check that the changes to the function interpid() are correct
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_ti = "InterpTurbine"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    # wind.input_shear == "Interpolation"
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_shear = "Interpolation"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    # yaw_method == "Constant"
    wind, sim, con, floris, floridyn = setup(settings_file)
    con.yaw = "Constant"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    # yaw_method == "InterpTurbine"
    wind, sim, con, floris, floridyn = setup(settings_file)
    con.yaw = "InterpTurbine"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)

    wind, sim, con, floris, floridyn = setup(settings_file)
    t = 0.0
    wind.input_dir = "RW_with_Mean"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)
    phi = getDataDir(set, wind, wf, t)

    wind, sim, con, floris, floridyn = setup(settings_file)
    t = 0.0
    wind.input_vel = "RW_with_Mean"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)
    # Create tmpM matrix with appropriate dimensions (nT × 3)
    # Column 1: velocity reduction factors, Column 2: added TI, Column 3: effective wind speed
    tmpM = [
        0.95 0.02 7.8;  # Turbine 1
        0.92 0.03 7.5;  # Turbine 2  
        0.90 0.04 7.3;  # Turbine 3
        0.94 0.02 7.7;  # Turbine 4
        0.88 0.05 7.2;  # Turbine 5
        0.85 0.06 6.9;  # Turbine 6
        0.93 0.02 7.6;  # Turbine 7
        0.87 0.05 7.1;  # Turbine 8
        0.82 0.07 6.7   # Turbine 9
    ]
    # Note: getDataVel call skipped for RW_with_Mean when wind.vel is Nothing
    # This combination requires proper wind.vel setup which is not needed for this test
    # u, wind = getDataVel(set, wind, wf, t, tmpM, floris)
end

@testset "readCovMatrix                                         " begin
    # Test scalar input (single value)
    @testset "Scalar input" begin
        nT = 3
        scalar_val = 0.5
        cov_mat, chol_factor = FLORIDyn.readCovMatrix([scalar_val], nT, "Test")
        
        @test size(cov_mat) == (nT, nT)
        @test cov_mat == scalar_val * I(nT)
        @test chol_factor == sqrt(scalar_val) * I(nT)
        @test all(diag(cov_mat) .== scalar_val)
        @test sum(cov_mat - Diagonal(cov_mat)) == 0  # Off-diagonal should be zero
    end
    
    # Test 1D vector input
    @testset "1D Vector input" begin
        nT = 3
        vec_vals = [0.1, 0.2, 0.3]
        cov_mat, chol_factor = FLORIDyn.readCovMatrix(vec_vals, nT, "Test")
        
        @test size(cov_mat) == (nT, nT)
        @test cov_mat == diagm(vec_vals)
        @test chol_factor == diagm(sqrt.(vec_vals))
        @test diag(cov_mat) == vec_vals
        @test sum(cov_mat - Diagonal(cov_mat)) == 0  # Off-diagonal should be zero
    end
    
    # Test row vector input (1×n matrix)
    @testset "Row vector input" begin
        nT = 3
        row_vec = [0.1 0.2 0.3]  # 1×3 matrix
        cov_mat, chol_factor = FLORIDyn.readCovMatrix(row_vec, nT, "Test")
        
        @test size(cov_mat) == (nT, nT)
        @test cov_mat == diagm(vec(row_vec))
        @test chol_factor == diagm(sqrt.(vec(row_vec)))
        @test diag(cov_mat) == vec(row_vec)
        @test sum(cov_mat - Diagonal(cov_mat)) == 0  # Off-diagonal should be zero
    end
    
    # Test column vector input (n×1 matrix)
    @testset "Column vector input" begin
        nT = 3
        col_vec = [0.1; 0.2; 0.3]  # 3×1 matrix
        cov_mat, chol_factor = FLORIDyn.readCovMatrix(col_vec, nT, "Test")
        
        @test size(cov_mat) == (nT, nT)
        @test cov_mat == diagm(vec(col_vec))
        @test chol_factor == diagm(sqrt.(vec(col_vec)))
        @test diag(cov_mat) == vec(col_vec)
        @test sum(cov_mat - Diagonal(cov_mat)) == 0  # Off-diagonal should be zero
    end
    
    # Test full matrix input (n×n matrix)
    @testset "Full matrix input" begin
        nT = 2
        full_mat = [0.1 0.05; 0.05 0.2]  # Symmetric positive definite matrix
        cov_mat, chol_factor = FLORIDyn.readCovMatrix(full_mat, nT, "Test")
        
        @test size(cov_mat) == (nT, nT)
        @test cov_mat == full_mat
        @test size(chol_factor) == (nT, nT)
        # Verify Cholesky decomposition: L * L' ≈ original matrix
        @test chol_factor * chol_factor' ≈ full_mat rtol=1e-10
        @test istril(chol_factor)  # Lower triangular
    end
    
    # Test error conditions
    @testset "Error conditions" begin
        nT = 3
        
        # Wrong vector size
        wrong_size_vec = [0.1, 0.2]  # Should be length 3
        @test_throws ErrorException FLORIDyn.readCovMatrix(wrong_size_vec, nT, "Test")
        
        # Wrong matrix dimensions
        wrong_mat = [0.1 0.2; 0.3 0.4; 0.5 0.6]  # 3×2 instead of 3×3
        @test_throws ErrorException FLORIDyn.readCovMatrix(wrong_mat, nT, "Test")
        
        # Non-square matrix with wrong first dimension
        wrong_mat2 = [0.1 0.2 0.3; 0.4 0.5 0.6]  # 2×3 instead of 3×3
        @test_throws ErrorException FLORIDyn.readCovMatrix(wrong_mat2, nT, "Test")
    end
    
    # Test edge cases
    @testset "Edge cases" begin
        # Single turbine
        nT = 1
        single_val = 0.25
        cov_mat, chol_factor = FLORIDyn.readCovMatrix([single_val], nT, "Test")
        @test size(cov_mat) == (1, 1)
        @test cov_mat[1,1] == single_val
        @test chol_factor[1,1] == sqrt(single_val)
        
        # Small positive variance (edge case that works)
        nT = 2
        small_val = 1e-6
        small_mat = [small_val 0.0; 0.0 small_val]
        cov_mat, chol_factor = FLORIDyn.readCovMatrix(small_mat, nT, "Test")
        @test cov_mat == small_mat
        @test chol_factor ≈ [sqrt(small_val) 0.0; 0.0 sqrt(small_val)] rtol=1e-10
    end
    
    # Test real-world scenarios
    @testset "Real-world scenarios" begin
        # Typical wind velocity covariance
        nT = 4
        wind_vel_cov = [0.5]  # Same variance for all turbines
        cov_mat, chol_factor = FLORIDyn.readCovMatrix(wind_vel_cov, nT, "WindVel")
        @test all(diag(cov_mat) .== 0.5)
        @test all(chol_factor[diagind(chol_factor)] .== sqrt(0.5))
        
        # Individual turbine variances
        individual_vars = [0.4, 0.6, 0.5, 0.7]
        cov_mat, chol_factor = FLORIDyn.readCovMatrix(individual_vars, nT, "WindVel")
        @test diag(cov_mat) == individual_vars
        @test diag(chol_factor) == sqrt.(individual_vars)
        
        # Correlated turbines (realistic covariance matrix)
        corr_mat = [0.5 0.1 0.05 0.02;
                    0.1 0.6 0.08 0.03;
                    0.05 0.08 0.4 0.04;
                    0.02 0.03 0.04 0.7]
        cov_mat, chol_factor = FLORIDyn.readCovMatrix(corr_mat, nT, "WindVel")
        @test cov_mat == corr_mat
        @test chol_factor * chol_factor' ≈ corr_mat rtol=1e-10
    end
end
nothing