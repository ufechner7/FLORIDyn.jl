using FLORIDyn
using Test

@testset "FLORIDyn.jl" begin
    WindDir = 270
    iT = [1, 2, 3]
    phi = getWindDirT(WindDir, iT, nothing)
    for ph in phi
        @test ph ≈ 270.0
    end
    
    # Suppose WindDir is a matrix where each row is [time, phi_T0, phi_T1, ...]
    WindDir = [
        0.0  10.0  20.0
        1.0  12.0  22.0
        2.0  14.0  24.0
    ]
    phi = getWindDirT_EnKF(WindDir, 1, 0.5)
    @test phi ≈ 11.0
end
