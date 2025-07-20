# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test

# Define a mock T object
mutable struct MockT
    States_WF::Matrix{Float64}
    StartI::Int
end

@testset "correctDir! updates turbine direction correctly" begin
    # Setup
    dir_strategy = Direction_All()
    mock_t = MockT(zeros(3, 4), 2)
    wind = :mockwind  # dummy input
    sim_time = :simtime  # dummy input

    # Call the function
    correctDir!(dir_strategy, mock_t, wind, sim_time)

    # Validate
    @test all(mock_t.States_WF[:, 2] .== 999.9)
    @test mock_t.States_WF[2, 4] == 999.9
end

@testset "correctDir! without 4th state column" begin
    dir_strategy = Direction_All()
    mock_t = MockT(zeros(3, 3), 1)
    wind = :mockwind
    sim_time = :simtime

    correctDir!(dir_strategy, mock_t, wind, sim_time)

    @test all(mock_t.States_WF[:, 2] .== 999.9)
    # Should not throw an error or modify column 4
    @test size(mock_t.States_WF, 2) == 3
end
