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

T, wind, sim, con, paramFLORIS = prepareSimulation(wind, con, paramFLORIDyn, paramFLORIS, turbProp, sim)

@testset "prepare_simulation                                      " begin
    @test wind.vel == 8.2
    @test wind.dir ≈ [
                        0.0    255.0
                    20600.0  255.0
                    20900.0  195.0
                    21200.0  195.0
                    ]
    @test wind.ti == 0.0620
end
nothing