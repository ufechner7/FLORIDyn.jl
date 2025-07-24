using FLORIDyn, Test

@testset verbose=true "settings" begin
    @testset "getTurbineData" begin
        # Test single known turbine types
        @testset "Single turbines" begin
            for (name, expected_pos, expected_D) in [
                ("DTU 10MW", [0.0, 0.0, 119.0], 178.4),
                ("DTU 5MW", [0.0, 0.0, 119.0], 178.4),
                ("Senvion 6.2M", [0.0, 0.0, 123.0], 126.0),
                ("V116", [0.0, 0.0, 84.0], 116.0),
                ("V117", [0.0, 0.0, 84.0], 117.0),
                ("V162", [0.0, 0.0, 119.0], 162.0),
                ("GE Haliade X", [0.0, 0.0, 150.0], 220.0)
            ]
                data = getTurbineData([name])
                @test size(data.NacPos) == (1,3)
                @test data.NacPos[1, :] ≈ expected_pos
                @test data.D[1] ≈ expected_D
            end
        end

        # Test multiple turbines
        @testset "Multiple turbines" begin
            names = ["DTU 10MW", "V117", "GE Haliade X"]
            data = getTurbineData(names)
            @test size(data.NacPos) == (3,3)
            @test data.NacPos[1, :] ≈ [0.0, 0.0, 119.0]
            @test data.D[1] ≈ 178.4
            @test data.NacPos[2, :] ≈ [0.0, 0.0, 84.0]
            @test data.D[2] ≈ 117.0
            @test data.NacPos[3, :] ≈ [0.0, 0.0, 150.0]
            @test data.D[3] ≈ 220.0
        end

        # # Test unknown turbine raises error
        # @testset "Unknown turbine error" begin
        #     err = @test_throws ErrorException getTurbineData(["Unknown Turbine"])
        #     println(err.msg)
        #     @test occursin("not known or misspelled", err.msg)
        # end

        # Test empty input returns empty arrays
        @testset "Empty input" begin
            data = getTurbineData(String[])
            @test size(data.NacPos) == (0, 3)
            @test length(data.D) == 0
        end
    end
end
nothing