using FLORIDyn
using Test
using LinearAlgebra
using Random

@testset "windshear" begin
    shear_mode = Shear_Interpolation()
    WindShear = [10 0.8; 20 0.9; 30 1.0]  # Example data
    z = [5, 15, 25, 35]                   # Heights to interpolate at
    shear = getWindShearT(shear_mode, WindShear, z)
    @test shear[1] ≈ 0.8
    @test shear[2] ≈ 0.85
    @test shear[3] ≈ 0.95
    @test shear[4] ≈ 1.0
end
nothing