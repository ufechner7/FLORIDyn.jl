@testset verbose=true "floris" begin
    @testset "discretizeRotor" begin

        # Test 1: Check output sizes
        nRP = 27  # example input
        RPs, w = discretizeRotor(nRP)
        N1 = 3
        n = round(Int, sqrt(nRP / N1))
        expected_nC = N1 * n^2

        @test size(RPs) == (expected_nC, 3)
        @test length(w) == expected_nC

        # Test 2: Check that RPs values in the 2nd and 3rd columns are in [-0.5, 0.5]
        @test all(abs.(RPs[:,2]) .<= 0.5)
        @test all(abs.(RPs[:,3]) .<= 0.5)

        # Test 3: Check that first column of RPs is all zeros (matches MATLAB code)
        @test all(RPs[:,1] .== 0.0)

        # Test 4: Check that weights sum to approximately 1
        @test isapprox(sum(w), 1.0; atol=1e-14)

        # Test 5: For different input, check outputs are plausible
        nRP2 = 48
        RPs2, w2 = discretizeRotor(nRP2)
        n2 = round(Int, sqrt(nRP2 / N1))
        @test size(RPs2) == (N1*n2^2, 3)
        @test all(abs.(RPs2[:,2]) .<= 0.5)
        @test all(abs.(RPs2[:,3]) .<= 0.5)
        @test isapprox(sum(w2), 1.0; atol=1e-14)
    end
    @testset "CalcCt" begin
        # Test typical values
        @test isapprox(CalcCt(0.0, nothing), 0.0)
        @test isapprox(CalcCt(0.25, nothing), 4 * 0.25 * 0.75)
        @test isapprox(CalcCt(0.5, nothing), 4 * 0.5 * 0.5)
        
        # Test edge cases
        @test isapprox(CalcCt(1.0, nothing), 0.0)
        @test isapprox(CalcCt(-0.1, nothing), 4 * -0.1 * (1 + 0.1))  # input out of typical range
        
        # Test vectorized input
        a_values = [0.0, 0.1, 0.2, 0.3]
        expected = 4 .* a_values .* (1 .- a_values)
        @test all(isapprox.(CalcCt.(a_values, Ref(nothing)), expected))
    end
end