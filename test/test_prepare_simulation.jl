# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test



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
    wind, sim, con, floris, floridyn = setup(settings_file)
    wind.input_dir = "RW_with_Mean"
    set = Settings(wind, sim, con)
    turbine_properties         = turbineArrayProperties(settings_file)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, turbine_properties, sim)
end
nothing