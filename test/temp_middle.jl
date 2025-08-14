@testset verbose=true "Pretty Print Functions                                  " begin
    @testset "turbines(T::Dict)" begin
        @testset "basic functionality with States_T (capital T)" begin
            # Test with the standard key format used in MAT files
            T_dict = Dict(
                "States_T" => [1.0 2.0 3.0; 4.0 5.0 6.0],  # (nT*nOP) × n_state_variables = (2*1) × 3
                "Names_T" => ["a", "yaw", "TI"],
                "nT" => 2,  # 2 turbines  
                "nOP" => 1  # 1 operating point each
            )
            
            df = turbines(T_dict)
            
            # Test output structure - now with OP and Turbine columns plus the data columns
            @test isa(df, DataFrame)
            @test size(df) == (2, 5)  # 2 turbines × 1 OP = 2 rows, 5 columns (OP + Turbine + 3 data)
            @test names(df) == ["OP", "Turbine", "a", "yaw", "TI"]
            
            # Test identifier columns
            @test df.OP == [1, 1]  # 1 OP for each turbine
            @test df.Turbine == [1, 2]  # Turbine numbers
            
            # Test data values - data rows are: turbine 1 OP 1, turbine 2 OP 1
            @test df.a == [1.0, 4.0]
            @test df.yaw == [2.0, 5.0]
            @test df.TI == [3.0, 6.0]
            
            # Test column types
            @test eltype(df.OP) == Int64
            @test eltype(df.Turbine) == Int64
            @test eltype(df.a) == Float64
            @test eltype(df.yaw) == Float64
            @test eltype(df.TI) == Float64
        end
        
        @testset "basic functionality with States_t (lowercase t)" begin
            # Test error when using lowercase variant (should fail)
            T_dict = Dict(
                "States_t" => [10.0 20.0; 30.0 40.0],  # This key should fail
                "Names_T" => ["T1", "T2"],
                "nT" => 2,
                "nOP" => 1
            )
            
            # This should throw an error since lowercase is not supported
            @test_throws ArgumentError turbines(T_dict)
        end
        
        @testset "symbol keys support" begin
            # Test with symbol keys instead of string keys
            T_dict = Dict(
                :States_T => [1.5 2.5; 3.5 4.5],  # 2 turbines × 1 OP each, 2 states
                :Names_T => ["A", "B"],
                :nT => 2,
                :nOP => 1
            )
            
            df = turbines(T_dict)
            
            @test isa(df, DataFrame)
            @test size(df) == (2, 4)  # 2 turbines × 1 OP = 2 rows, 4 columns
            @test names(df) == ["OP", "Turbine", "A", "B"]
            @test df.OP == [1, 1]
            @test df.Turbine == [1, 2]
            @test df.A == [1.5, 3.5]
            @test df.B == [2.5, 4.5]
        end
        
        @testset "mixed case key support" begin
            # Test mixed string/symbol keys
            T_dict = Dict(
                "States_T" => [0.1 0.2; 0.3 0.4],
                :Names_T => ["Mixed1", "Mixed2"],
                "nT" => 2,
                :nOP => 1
            )
            
            df_mixed = turbines(T_dict)
            
            @test isa(df_mixed, DataFrame)
            @test size(df_mixed) == (2, 4)  # 2 turbines × 1 OP = 2 rows, 4 columns
            @test names(df_mixed) == ["OP", "Turbine", "Mixed1", "Mixed2"]
        end
        
        @testset "single turbine" begin
            # Test with only one turbine
            T_dict = Dict(
                "States_T" => reshape([5.5, 6.6, 7.7], 1, 3),  # 1 turbine × 1 OP, 3 states
                "Names_T" => ["state1", "state2", "state3"],
                "nT" => 1,
                "nOP" => 1
            )
            
            df = turbines(T_dict)
            
            @test isa(df, DataFrame)
            @test size(df) == (1, 5)  # 1 turbine × 1 OP = 1 row, 5 columns  
            @test names(df) == ["OP", "Turbine", "state1", "state2", "state3"]
            @test df.OP == [1]
            @test df.Turbine == [1] 
            @test df.state1 == [5.5]
            @test df.state2 == [6.6]
            @test df.state3 == [7.7]
        end
        
        @testset "large dataset" begin
            # Test with larger dataset (multiple OPs per turbine)
            nT, nOP = 3, 4
            states_data = rand(nT * nOP, 2)  # (3×4) × 2 = 12 × 2 matrix
            T_dict = Dict(
                "States_T" => states_data,
                "Names_T" => ["param1", "param2"],
                "nT" => nT,
                "nOP" => nOP  
            )
            
            df = turbines(T_dict)
            
            @test isa(df, DataFrame)
            @test size(df) == (12, 4)  # 3 turbines × 4 OPs = 12 rows, 4 columns
            @test names(df) == ["OP", "Turbine", "param1", "param2"]
            
            # Check OP pattern: 1,2,3,4,1,2,3,4,1,2,3,4 (4 OPs repeated for 3 turbines)
            @test df.OP == repeat(1:nOP, nT)
            # Check Turbine pattern: 1,1,1,1,2,2,2,2,3,3,3,3 (each turbine repeated 4 times)  
            @test df.Turbine == repeat(1:nT, inner=nOP)
        end
        
        @testset "different data types" begin
            # Test with different numeric types
            T_dict = Dict(
                "States_T" => Int32[10 20; 30 40],  # Integer input
                "Names_T" => ["IntVal1", "IntVal2"],
                "nT" => 2,
                "nOP" => 1
            )
            
            df = turbines(T_dict)
            
            @test isa(df, DataFrame)
            @test size(df) == (2, 4)
            @test names(df) == ["OP", "Turbine", "IntVal1", "IntVal2"] 
            @test df.IntVal1 == [10, 30]
            @test df.IntVal2 == [20, 40]
        end
        
        @testset "real MAT file format" begin
            # Test with format similar to actual MAT file data
            nT, nOP = 3, 2
            T_dict = Dict(
                "States_T" => [
                    0.33 -0.1 0.0;      # Turbine 1 OP 1
                    0.33 -0.099 0.001;  # Turbine 1 OP 2  
                    0.33 -0.1 0.0;      # Turbine 2 OP 1
                    0.33 -0.099 0.001;  # Turbine 2 OP 2
                    0.33 -0.1 0.0;      # Turbine 3 OP 1
                    0.33 -0.099 0.001   # Turbine 3 OP 2
                ],  # 6 rows = nT*nOP, 3 cols = n_state_variables
                "Names_T" => ["a", "yaw", "TI"],  # Actual names from MAT files
                "nT" => nT,
                "nOP" => nOP
            )
            
            df = turbines(T_dict)
            
            @test isa(df, DataFrame)
            @test size(df) == (6, 5)  # 3 turbines × 2 OPs = 6 rows, 5 columns
            @test names(df) == ["OP", "Turbine", "a", "yaw", "TI"]
            @test df.OP == [1, 2, 1, 2, 1, 2]  # Pattern: OP 1-2 for each turbine
            @test df.Turbine == [1, 1, 2, 2, 3, 3]  # Pattern: each turbine repeated nOP times
            @test all(df.a .== 0.33)  # All axial induction values same
        end
        
        @testset "error handling - missing keys" begin
            # Test error when States key is missing
            @test_throws ArgumentError turbines(Dict("Names_T" => ["T1"], "nT" => 1, "nOP" => 1))
            @test_throws ArgumentError turbines(Dict("Random_Key" => [1.0 2.0], "Names_T" => ["T1", "T2"], "nT" => 2, "nOP" => 1))
            
            # Test error when Names key is missing  
            @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "nT" => 2, "nOP" => 1))
            @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "Random_Key" => ["T1", "T2"], "nT" => 2, "nOP" => 1))
            
            # Test error when nT key is missing
            @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "Names_T" => ["T1", "T2"], "nOP" => 1))
            
            # Test error when nOP key is missing
            @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "Names_T" => ["T1", "T2"], "nT" => 2))
            
            # Test error when using lowercase keys (not supported)
            @test_throws ArgumentError turbines(Dict("States_t" => [1.0 2.0; 3.0 4.0], "Names_T" => ["T1", "T2"], "nT" => 2, "nOP" => 1))
            @test_throws ArgumentError turbines(Dict("States_T" => [1.0 2.0; 3.0 4.0], "Names_t" => ["T1", "T2"], "nT" => 2, "nOP" => 1))
        end
        
        @testset "error handling - empty data" begin
            # Test error with empty states matrix
            @test_throws ArgumentError turbines(Dict(
                "States_T" => Matrix{Float64}(undef, 0, 0),
                "Names_T" => String[],
                "nT" => 0,
                "nOP" => 0
            ))
            
            # Test error with empty names vector
            @test_throws ArgumentError turbines(Dict(
                "States_T" => [1.0 2.0; 3.0 4.0],
                "Names_T" => String[],
                "nT" => 2,
                "nOP" => 1
            ))
        end
        
        @testset "error handling - dimension mismatch" begin
            # Test error when number of columns doesn't match number of names
            @test_throws DimensionMismatch turbines(Dict(
                "States_T" => [1.0 2.0 3.0; 4.0 5.0 6.0],  # 3 state variables
                "Names_T" => ["T1", "T2"],  # Only 2 names
                "nT" => 2,
                "nOP" => 1
            ))
            
            @test_throws DimensionMismatch turbines(Dict(
                "States_T" => [1.0 2.0; 3.0 4.0],  # 2 state variables
                "Names_T" => ["T1", "T2", "T3"],  # 3 names
                "nT" => 2,
                "nOP" => 1
            ))
            
            # Test error when States_T dimensions don't match nT*nOP
            @test_throws DimensionMismatch turbines(Dict(
                "States_T" => [1.0 2.0; 3.0 4.0],  # 2 rows 
                "Names_T" => ["T1", "T2"],  
                "nT" => 3,  # 3 turbines
                "nOP" => 2  # 2 OPs -> should be 6 rows, but only have 2
            ))
        end
        
        @testset "error handling - invalid dictionary" begin
            # Test with empty dictionary
            @test_throws ArgumentError turbines(Dict{String, Any}())
        end
        
        @testset "integration with actual MAT file data" begin
            # Test that the function works with the format from test_case_01.jl
            # This simulates the actual MAT file structure
            using MAT  # Ensure MAT is available for this test
            
            # Create test data similar to what would come from matread()
            nT, nOP = 2, 3
            T_ref_simulation = Dict(
                "States_T" => [
                    0.33  -0.09999  0.0;   # Turbine 1, OP 1
                    0.33  -0.1      0.0;   # Turbine 1, OP 2
                    0.33  -0.099    0.05;  # Turbine 1, OP 3
                    0.33  -0.10001  0.0;   # Turbine 2, OP 1
                    0.33  -0.1      0.0;   # Turbine 2, OP 2
                    0.33  -0.098    0.03   # Turbine 2, OP 3
                ],
                "Names_T" => ["a", "yaw", "TI"],  # State variable names
                "nT" => nT,
                "nOP" => nOP
            )
            
            df = turbines(T_ref_simulation)
            
            @test isa(df, DataFrame)
            @test size(df) == (6, 5)  # 2 turbines × 3 OPs = 6 rows, 5 columns
            @test names(df) == ["OP", "Turbine", "a", "yaw", "TI"]
            
            # Verify the structure matches expected patterns
            @test df.OP == [1, 2, 3, 1, 2, 3]  # 1,2,3 for turbine 1, then 1,2,3 for turbine 2
            @test df.Turbine == [1, 1, 1, 2, 2, 2]  # Each turbine repeated for all OPs
            @test all(df.a .== 0.33)  # Axial induction should be constant
            @test df.TI[3] ≈ 0.05  # OP 3 for turbine 1
            @test df.TI[6] ≈ 0.03  # OP 3 for turbine 2
        end
        
        @testset "function integration" begin
            # Test that the function is properly exported
            @test :turbines in names(FLORIDyn)
            
            # Test it works with the same data format as turbines(wf::WindFarm)
            T_dict = Dict(
                "States_T" => [0.33 -0.1 0.06; 0.33 -0.1 0.06],  # 2 turbines × 1 OP, 3 states
                "Names_T" => ["a", "yaw", "TI"],
                "nT" => 2,
                "nOP" => 1
            )
            
            df = turbines(T_dict)
            
            # Should have same columns as WindFarm version
            @test "OP" in names(df)
            @test "Turbine" in names(df) 
            @test "a" in names(df)
            @test "yaw" in names(df)
            @test "TI" in names(df)
        end
        
        @testset "edge cases" begin
            # Test with minimal valid data (1 turbine, 1 OP, 1 state)
            T_dict = Dict(
                "States_T" => reshape([42.0], 1, 1),
                "Names_T" => ["single_state"],
                "nT" => 1,
                "nOP" => 1
            )
            
            df = turbines(T_dict)
            @test isa(df, DataFrame)
            @test size(df) == (1, 3)  # 1 row, 3 columns (OP + Turbine + data)
            @test df.OP == [1]
            @test df.Turbine == [1]
            @test df.single_state == [42.0]
        end
        
        @testset "performance with large datasets" begin
            # Test that the function performs reasonably with large data
            nT, nOP = 10, 100
            n_states = 5
            
            large_states = rand(Float64, nT * nOP, n_states)  # 1000 × 5 matrix
            large_names = ["State_$i" for i in 1:n_states]
            
            T_large = Dict(
                "States_T" => large_states,
                "Names_T" => large_names,
                "nT" => nT,
                "nOP" => nOP
            )
            
            # Time the operation (should complete in reasonable time)
            time_start = time()
            df_large = turbines(T_large)
            time_elapsed = time() - time_start
            
            @test isa(df_large, DataFrame)
            @test size(df_large) == (1000, 7)  # 1000 rows × (2 + 5) columns
            @test time_elapsed < 5.0  # Should complete in less than 5 seconds
            
            # Verify data integrity for large dataset
            @test df_large.OP == repeat(1:nOP, nT)  # Pattern repeats correctly
            @test df_large.Turbine == repeat(1:nT, inner=nOP)  # Turbine pattern is correct
        end
    end
end
