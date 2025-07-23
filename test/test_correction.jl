# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test

correction = FLORIDyn.WindCorrection("None", "All", "None")
pertubation = FLORIDyn.WindPerturbation(0.0, 0.2, 0.0, 0.5, 0.0, 0.005)
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
    pertubation,
    8.2,
    dir_array,
    0.062,
    shear
)
@testset verbose=true "correction" begin
    @testset "correctDir!" begin
        # Setup
        set = Settings(Velocity_Constant(), Direction_Interpolation(), TI_Constant(), Shear_PowerLaw(), Direction_All(), 
                    Velocity_None(), TI_None(), IterateOPs_basic(), Yaw_SOWFA())
        T = WindFarm(
            States_WF = zeros(3, 4),
            StartI = [2 2; 3 3],
            nT = 3,
            D = [126.0, 126.0, 126.0],
            nOP = 4,
            intOPs = [[1 0.5 2 0.5; 2 0.5 3 0.5]]
        )
         sim_time = 20000  # dummy input

        # Call the function
        correctDir!(set.dir_mode, set, T, wind, sim_time)

        # Validate
        @test all(T.States_WF[:, 2] .== 255.0)
        @test T.States_WF[2, 4] == 255.0
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