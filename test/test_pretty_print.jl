# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Unit tests for pretty_print.jl functionality

if !isdefined(Main, :Test)
    using Test
end

if !isdefined(Main, :FLORIDyn)
    using FLORIDyn
end

if !isdefined(Main, :DataFrames)
    using DataFrames
end

@testset verbose=true "Pretty Print Tests" begin

    # Helper function to create a test WindFarm
    function create_test_windfarm(;nT=3, nOP=24, nWF=2, nStates=4, nTimesteps=12)
        WindFarm(
            nT = nT,
            nOP = nOP,
            posBase = rand(2, nT),
            D = fill(126.0, nT),
            Names_T = ["T$i" for i in 1:nT],
            Names_WF = ["WF_$i" for i in 1:nWF],
            Names_OP = ["OP_$i" for i in 1:3],  # Always 3 OP variables for simplicity
            States_WF = rand(nTimesteps, nWF),  # nTimesteps rows × nWF columns  
            States_T = rand(nStates, nT),       # nStates rows × nT columns
            States_OP = rand(nOP, 3),          # nOP rows × 3 columns (always 3 OP vars)
            dep = [Int[] for _ in 1:nT]
        )
    end

    @testset "turbines() function" begin
        @testset "Basic functionality" begin
            wf = create_test_windfarm()
            df = FLORIDyn.turbines(wf)
            
            @test df isa DataFrame
            @test size(df, 1) == 4  # 4 states per turbine
            @test size(df, 2) == 3  # 3 turbines
            @test names(df) == ["T1", "T2", "T3"]
            @test eltype(df.T1) == Float64
        end
        
        @testset "Single turbine" begin
            wf = create_test_windfarm(nT=1, nStates=5)
            df = FLORIDyn.turbines(wf)
            
            @test size(df) == (5, 1)
            @test names(df) == ["T1"]
        end
        
        @testset "Many turbines" begin
            wf = create_test_windfarm(nT=10, nStates=6)
            df = FLORIDyn.turbines(wf)
            
            @test size(df) == (6, 10)
            @test length(names(df)) == 10
            @test names(df)[1] == "T1"
            @test names(df)[10] == "T10"
        end
        
        @testset "Error handling" begin
            # Empty States_T
            @test_throws ArgumentError begin
                wf = WindFarm(
                    nT = 2, nOP = 10, posBase = rand(2, 2), D = [126.0, 126.0],
                    Names_T = String[], Names_WF = ["Power"], Names_OP = ["WindSpeed"],
                    States_WF = rand(5, 1), States_T = Matrix{Float64}(undef, 0, 0),
                    States_OP = rand(10, 1), dep = [Int[], Int[]]
                )
                FLORIDyn.turbines(wf)
            end
            
            # Dimension mismatch
            @test_throws DimensionMismatch begin
                wf = WindFarm(
                    nT = 2, nOP = 10, posBase = rand(2, 2), D = [126.0, 126.0],
                    Names_T = ["T1", "T2", "T3"], Names_WF = ["Power"], Names_OP = ["WindSpeed"],
                    States_WF = rand(5, 1), States_T = rand(3, 2), States_OP = rand(10, 1),
                    dep = [Int[], Int[]]
                )
                FLORIDyn.turbines(wf)
            end
        end
    end

    @testset "windfield() function" begin
        @testset "Basic functionality" begin
            wf = create_test_windfarm(nWF=3, nTimesteps=8)
            df = FLORIDyn.windfield(wf)
            
            @test df isa DataFrame
            @test size(df, 1) == 8   # 8 time steps
            @test size(df, 2) == 3   # 3 wind field variables
            @test names(df) == ["WF_1", "WF_2", "WF_3"]
            @test eltype(df.WF_1) == Float64
        end
        
        @testset "Single wind field variable" begin
            wf = create_test_windfarm(nWF=1, nTimesteps=10)
            df = FLORIDyn.windfield(wf)
            
            @test size(df) == (10, 1)
            @test names(df) == ["WF_1"]
        end
        
        @testset "Many timesteps" begin
            wf = create_test_windfarm(nWF=2, nTimesteps=100)
            df = FLORIDyn.windfield(wf)
            
            @test size(df) == (100, 2)
            @test names(df) == ["WF_1", "WF_2"]
        end
        
        @testset "Error handling" begin
            # Empty States_WF
            @test_throws ArgumentError begin
                wf = WindFarm(
                    nT = 2, nOP = 10, posBase = rand(2, 2), D = [126.0, 126.0],
                    Names_T = ["T1", "T2"], Names_WF = String[], Names_OP = ["WindSpeed"],
                    States_WF = Matrix{Float64}(undef, 0, 0), States_T = rand(3, 2),
                    States_OP = rand(10, 1), dep = [Int[], Int[]]
                )
                FLORIDyn.windfield(wf)
            end
            
            # Dimension mismatch
            @test_throws DimensionMismatch begin
                wf = WindFarm(
                    nT = 2, nOP = 10, posBase = rand(2, 2), D = [126.0, 126.0],
                    Names_T = ["T1", "T2"], Names_WF = ["Power", "Thrust", "Efficiency"],
                    Names_OP = ["WindSpeed"], States_WF = rand(5, 2), States_T = rand(3, 2),
                    States_OP = rand(10, 1), dep = [Int[], Int[]]
                )
                FLORIDyn.windfield(wf)
            end
        end
    end

    @testset "ops() function" begin
        @testset "Basic functionality" begin
            wf = create_test_windfarm(nOP=20)
            df = FLORIDyn.ops(wf)
            
            @test df isa DataFrame
            @test size(df, 1) == 20  # 20 operating points
            @test size(df, 2) == 3   # 3 operating point variables
            @test names(df) == ["OP_1", "OP_2", "OP_3"]
            @test eltype(df.OP_1) == Float64
        end
        
        @testset "Single operating point variable" begin
            wf = WindFarm(
                nT = 2, nOP = 15, posBase = rand(2, 2), D = [126.0, 126.0],
                Names_T = ["T1", "T2"], Names_WF = ["Power", "Thrust"],
                Names_OP = ["WindSpeed"], States_WF = rand(8, 2), States_T = rand(3, 2),
                States_OP = rand(15, 1), dep = [Int[], Int[]]
            )
            df = FLORIDyn.ops(wf)
            
            @test size(df) == (15, 1)
            @test names(df) == ["WindSpeed"]
        end
        
        @testset "Vector States_OP" begin
            wf = WindFarm(
                nT = 2, nOP = 10, posBase = rand(2, 2), D = [126.0, 126.0],
                Names_T = ["T1", "T2"], Names_WF = ["Power", "Thrust"],
                Names_OP = ["WindSpeed"], States_WF = rand(5, 2), States_T = rand(3, 2),
                States_OP = reshape(rand(10), 10, 1), dep = [Int[], Int[]]  # Make it a matrix
            )
            df = FLORIDyn.ops(wf)
            
            @test size(df) == (10, 1)
            @test names(df) == ["WindSpeed"]
        end
        
        @testset "Error handling" begin
            # Empty States_OP
            @test_throws ArgumentError begin
                wf = WindFarm(
                    nT = 2, nOP = 10, posBase = rand(2, 2), D = [126.0, 126.0],
                    Names_T = ["T1", "T2"], Names_WF = ["Power"], Names_OP = String[],
                    States_WF = rand(5, 1), States_T = rand(3, 2),
                    States_OP = Matrix{Float64}(undef, 0, 0), dep = [Int[], Int[]]
                )
                FLORIDyn.ops(wf)
            end
            
            # Dimension mismatch - matrix case
            @test_throws DimensionMismatch begin
                wf = WindFarm(
                    nT = 2, nOP = 10, posBase = rand(2, 2), D = [126.0, 126.0],
                    Names_T = ["T1", "T2"], Names_WF = ["Power"], 
                    Names_OP = ["U", "Dir", "TI"], States_WF = rand(5, 1), States_T = rand(3, 2),
                    States_OP = rand(10, 2), dep = [Int[], Int[]]  # 2 cols but 3 names
                )
                FLORIDyn.ops(wf)
            end
            
            # Dimension mismatch - vector case
            @test_throws DimensionMismatch begin
                wf = WindFarm(
                    nT = 2, nOP = 10, posBase = rand(2, 2), D = [126.0, 126.0],
                    Names_T = ["T1", "T2"], Names_WF = ["Power"], 
                    Names_OP = ["U", "Dir"], States_WF = rand(5, 1), States_T = rand(3, 2),
                    States_OP = reshape(rand(10), 10, 1), dep = [Int[], Int[]]  # 1 col but 2 names
                )
                FLORIDyn.ops(wf)
            end
        end
    end

    @testset "Property accessors (getproperty)" begin
        @testset "wf.turbines property" begin
            wf = create_test_windfarm()
            df_prop = wf.turbines
            df_func = FLORIDyn.turbines(wf)
            
            @test df_prop == df_func
            @test df_prop isa DataFrame
            @test size(df_prop) == (4, 3)
        end
        
        @testset "wf.windfield property" begin
            wf = create_test_windfarm()
            df_prop = wf.windfield
            df_func = FLORIDyn.windfield(wf)
            
            @test df_prop == df_func
            @test df_prop isa DataFrame
            @test size(df_prop) == (12, 2)
        end
        
        @testset "wf.ops property" begin
            wf = create_test_windfarm()
            df_prop = wf.ops
            df_func = FLORIDyn.ops(wf)
            
            @test df_prop == df_func
            @test df_prop isa DataFrame
            @test size(df_prop) == (24, 3)
        end
        
        @testset "Normal properties still work" begin
            wf = create_test_windfarm()
            
            @test wf.nT == 3
            @test wf.nOP == 24
            @test wf.Names_T == ["T1", "T2", "T3"]
            @test wf.Names_WF == ["WF_1", "WF_2"]
            @test wf.Names_OP == ["OP_1", "OP_2", "OP_3"]
            @test length(wf.D) == 3
        end
    end

    @testset "Pretty printing functions" begin
        @testset "WindFarm show method" begin
            wf = create_test_windfarm()
            
            # Test that show doesn't error
            io = IOBuffer()
            show(io, wf)
            output = String(take!(io))
            
            @test length(output) > 0
            @test occursin("WindFarm:", output)
            @test occursin("Turbines:", output)
            @test occursin("Operating Points:", output)
        end
        
        @testset "WindFarm summary method" begin
            wf = create_test_windfarm()
            summary_str = summary(wf)
            
            @test summary_str isa String
            @test occursin("WindFarm:", summary_str)
            @test occursin("3 turbines", summary_str)
            @test occursin("24 op. points", summary_str)
        end
        
        @testset "Layout estimation" begin
            # Grid layout
            wf = WindFarm(
                nT = 9, nOP = 10,
                posBase = [0.0 400.0 800.0 0.0 400.0 800.0 0.0 400.0 800.0;
                          0.0 0.0 0.0 400.0 400.0 400.0 800.0 800.0 800.0],
                D = fill(126.0, 9), Names_T = ["T$i" for i in 1:9],
                Names_WF = ["Power"], Names_OP = ["WindSpeed"],
                States_WF = rand(5, 1), States_T = rand(3, 9), States_OP = rand(10, 1),
                dep = [Int[] for _ in 1:9]
            )
            
            layout = FLORIDyn.estimate_layout(wf.posBase)
            @test layout isa String
            @test length(layout) > 0
            
            # Linear layout
            wf_linear = WindFarm(
                nT = 4, nOP = 10,
                posBase = [0.0 400.0 800.0 1200.0; 0.0 0.0 0.0 0.0],
                D = fill(126.0, 4), Names_T = ["T$i" for i in 1:4],
                Names_WF = ["Power"], Names_OP = ["WindSpeed"],
                States_WF = rand(5, 1), States_T = rand(3, 4), States_OP = rand(10, 1),
                dep = [Int[] for _ in 1:4]
            )
            
            layout_linear = FLORIDyn.estimate_layout(wf_linear.posBase)
            @test layout_linear isa String
            @test length(layout_linear) > 0
        end
    end

    @testset "Integration tests" begin
        @testset "All properties work together" begin
            wf = create_test_windfarm(nT=5, nOP=50, nWF=4, nStates=7, nTimesteps=20)
            
            # Test all three properties
            turb_df = wf.turbines
            wf_df = wf.windfield  
            ops_df = wf.ops
            
            @test size(turb_df) == (7, 5)
            @test size(wf_df) == (20, 4)
            @test size(ops_df) == (50, 3)
            
            # Test column names
            @test names(turb_df) == ["T1", "T2", "T3", "T4", "T5"]
            @test names(wf_df) == ["WF_1", "WF_2", "WF_3", "WF_4"]
            @test names(ops_df) == ["OP_1", "OP_2", "OP_3"]
        end
        
        @testset "Data access patterns" begin
            wf = create_test_windfarm()
            
            # Individual columns
            @test wf.turbines.T1 isa Vector{Float64}
            @test wf.windfield.WF_1 isa Vector{Float64}
            @test wf.ops.OP_1 isa Vector{Float64}
            
            # Specific indexing
            @test wf.turbines[1, "T2"] isa Float64
            @test wf.windfield[2, "WF_2"] isa Float64
            @test wf.ops[3, "OP_3"] isa Float64
            
            # Row access
            @test size(wf.turbines[1, :]) == (3,)  # 3 turbines
            @test size(wf.windfield[1, :]) == (2,) # 2 wf variables  
            @test size(wf.ops[1, :]) == (3,)       # 3 op variables
        end
        
        @testset "Consistency across calls" begin
            wf = create_test_windfarm()
            
            # Multiple calls should return same data
            df1 = wf.turbines
            df2 = wf.turbines
            @test df1 == df2
            
            df3 = wf.windfield
            df4 = wf.windfield
            @test df3 == df4
            
            df5 = wf.ops
            df6 = wf.ops
            @test df5 == df6
        end
    end
end
