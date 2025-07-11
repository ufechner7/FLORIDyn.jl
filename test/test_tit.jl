# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test

@testset "wind turbulence" begin
    WindTi = 0.1
    iT = [1, 2, 3]
    Ti = getWindTiT(TI_Constant(), WindTi, iT)
    # Ti will be [0.1, 0.1, 0.1]
    @test length(Ti) == 3
    @test Ti[1] ≈ 0.1
    @test Ti[2] ≈ 0.1
    @test Ti[3] ≈ 0.1
end
