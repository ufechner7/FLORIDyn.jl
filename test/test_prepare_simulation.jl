# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn, Test

settings_file = "data/2021_9T_Data.yaml"

# get the settings for the wind field, simulator and controller
wind, sim, con, floris, floridyn = setup(settings_file)

# create settings struct
set = Settings(wind, sim)

# % Load linked data
turbProp        = turbineArrayProperties(settings_file)
paramFLORIS     = floris
paramFLORIDyn   = floridyn

T, wind, sim, con, paramFLORIS = prepareSimulation(set, wind, con, paramFLORIDyn, paramFLORIS, turbProp, sim)

@testset "prepare_simulation                                      " begin
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

    @test T[:posBase] == expected_posBase
    @test T[:nT] == 9
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

    @test T[:posNac]    == expected

    @test T[:D]         == fill(178.4, 9)
    @test T[:Names_OP]  == ["x0", "y0", "z0", "x1", "y1", "z1"]
    @test T[:Names_T]   == ["a", "yaw", "TI"]
    @test T[:Names_WF]  == ["wind_vel", "wind_dir", "TI0", "OP_ori"]
    @test T[:StartI]    == [1 201 401 601 801 1001 1201 1401 1601]
    @test T[:nOP]       == 200
    @test T[:red_arr]   == ones(9,9)

    @test paramFLORIDyn.deltaUW == 10.0
    @test size(con.yaw_data) == (603, 10)
    @test sum(con.yaw_data) ≈ 1.37333255e7
    @test con.tanh_yaw == false

    @test size(T[:States_OP]) == (1800, 6)
    @test sum(T[:States_OP]) ≈ 1.8683e+07 rtol=1e-4

    @test size(T[:States_T]) == (1800, 3)
    @test sum(T[:States_T]) ≈ 702 rtol=1e-4

    @test size(T[:States_WF]) == (1800, 4)
    @test sum(T[:States_WF]) ≈ 9.3287e+05 rtol=1e-4

    @test sim.n_sim_steps == 301
end
nothing