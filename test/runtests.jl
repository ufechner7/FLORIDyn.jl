using FLORIDyn
using Test
using LinearAlgebra
using Random

Random.seed!(1234)

struct WindDirType
    Data::Float64
    CholSig::Matrix{Float64}
end

@testset "FLORIDyn.jl" begin
    dir_mode = Direction_Constant()
    WindDir = 270
    iT = [1, 2, 3]
    phi = getWindDirT(dir_mode, WindDir, iT, nothing)
    for ph in phi
        @test ph ≈ 270.0
    end

    dir_mode = Direction_Constant_wErrCov()
    WindDir = WindDirType(270.0, cholesky(Matrix{Float64}(I, 3, 3)).L)
    iT = [1, 2, 3]
    phi = getWindDirT(dir_mode, WindDir, iT)
    result = [269.6402710931765, 271.0872084924286, 269.5804103830612]
    for (i, ph) in pairs(phi)
        @test ph ≈ result[i]
    end

    dir_mode = Direction_EnKF_InterpTurbine()

    # Suppose WindDir is a matrix where each row is [time, phi_T0, phi_T1, ...]
    WindDir = [
        0.0  10.0  20.0
        1.0  12.0  22.0
        2.0  14.0  24.0
    ]
    phi = getWindDirT_EnKF(dir_mode, WindDir, 1, 0.5)
    @test phi ≈ 11.0
end
