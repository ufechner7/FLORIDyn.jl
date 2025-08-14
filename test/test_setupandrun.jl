# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause
 
using FLORIDyn, Test, Statistics, MAT, DataFrames

matlab_file   = "test/data/input_setUpTmpWFAndRun_2_steps.mat"
vars = matread(matlab_file)
wf_dict = vars["T"]

"""
    convert_wf_dict2windfarm(wf_dict::Dict{String, Any}) -> WindFarm

Convert MATLAB wf_dict to Julia WindFarm type.

This function handles the conversion of WindFarm data structures from MATLAB MAT files
to the Julia WindFarm type. It properly handles the nested cell arrays, mixed data types,
and different array orientations that result from MATLAB to Julia conversion.

# Arguments
- `wf_dict::Dict{String, Any}`: Dictionary from MAT file containing WindFarm data with keys:
  - `"nT"`: Number of turbines (scalar)
  - `"nOP"`: Number of operating points (scalar) 
  - `"States_WF"`: Wind field states matrix
  - `"States_OP"`: Operating point states matrix
  - `"States_T"`: Turbine states matrix
  - `"posBase"`: Base positions matrix (nT × 3)
  - `"posNac"`: Nacelle positions matrix (nT × 3) 
  - `"D"`: Turbine diameters (nT × 1)
  - `"StartI"`: Start indices (1 × nT)
  - `"intOPs"`: Interpolated operating points (nT × 1 cell array)
  - `"weight"`: Weights (nT × 1 cell array, may contain scalars or matrices)
  - `"dep"`: Dependencies (nT × 1 cell array, may contain scalars, vectors or matrices)
  - `"red_arr"`: Reduced array matrix
  - `"Names_T"`: Turbine state names (1 × n_states)
  - `"Names_WF"`: Wind field names (1 × n_wf)
  - `"Names_OP"`: Operating point names (1 × n_op)

# Returns
- `WindFarm`: Properly constructed WindFarm object with all fields correctly typed

# Examples
```julia
# Load MAT file and convert
vars = matread("wind_farm_data.mat")
wf_dict = vars["wf"]
wf = convert_wf_dict2windfarm(wf_dict)

# Access as DataFrames
turbine_data = wf.turbines
windfield_data = wf.windfield  
ops_data = wf.ops
```

# Notes
- Handles mixed scalar/array types in weight and dep fields
- Transposes position matrices to match Julia convention (3 × nT)
- Converts MATLAB row vectors to Julia column vectors for names
- Ensures all nested structures are properly typed for WindFarm constructor
"""
function convert_wf_dict2windfarm(wf_dict::Dict{String, Any})
    # Extract and convert the complex nested structures
    nT = Int(wf_dict["nT"])
    nOP = Int(wf_dict["nOP"])
    
    # Convert nested cell arrays
    intOPs_vec = Vector{Matrix{Float64}}()
    weight_vec = Vector{Vector{Float64}}()
    dep_vec = Vector{Vector{Int}}()
    
    for i in 1:nT
        push!(intOPs_vec, wf_dict["intOPs"][i])
        
        # Convert weight (may be empty matrix, scalar, or vector)
        w = wf_dict["weight"][i]
        if isempty(w)
            push!(weight_vec, Float64[])
        elseif isa(w, Number)  # Handle scalar case
            push!(weight_vec, [Float64(w)])
        else
            push!(weight_vec, vec(Float64.(w)))
        end
        
        # Convert dep (may be empty matrix, scalar, or vector)
        d = wf_dict["dep"][i]
        if isempty(d)
            push!(dep_vec, Int[])
        elseif isa(d, Number)  # Handle scalar case
            push!(dep_vec, [Int(d)])
        else
            push!(dep_vec, vec(Int.(d)))
        end
    end
    
    return WindFarm(
        nT = nT,
        nOP = nOP,
        States_WF = wf_dict["States_WF"],
        States_OP = wf_dict["States_OP"], 
        States_T = wf_dict["States_T"],
        posBase = wf_dict["posBase"]',  # Transpose to get (2 or 3) × nT matrix
        posNac = wf_dict["posNac"]',   # Transpose to get (2 or 3) × nT matrix  
        D = vec(wf_dict["D"]),         # Ensure it's a vector
        StartI = Int.(wf_dict["StartI"]),
        intOPs = intOPs_vec,
        Weight = weight_vec,
        dep = dep_vec,
        red_arr = wf_dict["red_arr"],
        Names_T = String.(wf_dict["Names_T"][1, :]),      # Get row vector and convert to strings
        Names_WF = String.(wf_dict["Names_WF"][1, :]),    # Get row vector and convert to strings
        Names_OP = String.(wf_dict["Names_OP"][1, :])     # Get row vector and convert to strings
    )
end

function structs_equal(a::T, b::T; prn=true) where T
    result = true
    fields = fieldnames(T)
    for f in fields
        val_a = getfield(a, f)
        val_b = getfield(b, f)
        if val_a != val_b
            prn && println("Field $(f): a = $(val_a), b = $(val_b)")
            result = false
        end
    end
    return result
end

 @testset "setUpTmpWFAndRun_basic" begin
    settings_file = "data/2021_9T_Data.yaml"
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    # create settings struct
    set = Settings(wind, sim, con)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    wf.dep = findTurbineGroups(wf, floridyn)
    wf.intOPs = interpolateOPs(wf)
    wf_old = deepcopy(wf)
    M, wf = setUpTmpWFAndRun(set, wf, floris, wind)
    @test ! structs_equal(wf_old, wf; prn=false)
end

@testset "wf_dict_to_windfarm_conversion" begin
    # Test the conversion function with MAT file data
    wf_converted = convert_wf_dict2windfarm(wf_dict)
    
    @test wf_converted isa WindFarm
    @test wf_converted.nT == 9
    @test wf_converted.nOP == 200
    @test wf_converted.Names_T == ["a", "yaw", "TI"]
    @test wf_converted.Names_WF == ["wind_vel", "wind_dir", "TI0", "OP_ori"] 
    @test wf_converted.Names_OP == ["x0", "y0", "z0", "x1", "y1", "z1"]
    @test size(wf_converted.posBase) == (3, 9)
    @test length(wf_converted.D) == 9
    @test length(wf_converted.intOPs) == 9
    @test length(wf_converted.Weight) == 9
    @test length(wf_converted.dep) == 9
    
    # Test DataFrame property access
    turbines_df = wf_converted.turbines
    @test turbines_df isa DataFrame
    @test size(turbines_df) == (1800, 5)  # nOP×nT = 200×9 = 1800 rows, 5 columns (OP + Turbine + 3 states)
    @test names(turbines_df) == ["OP", "Turbine", "a", "yaw", "TI"]
    
    wf_df = wf_converted.windfield
    @test wf_df isa DataFrame
    @test size(wf_df) == (1800, 4)  # 1800 timesteps, 4 wind field variables
    
    ops_df = wf_converted.ops
    @test ops_df isa DataFrame
    @test size(ops_df) == (1800, 6)  # 1800 states, 6 OP variables
end
 @testset "setUpTmpWFAndRun_vs_matlab" begin
    settings_file = "data/2021_9T_Data.yaml"
    # get the settings for the wind field, simulator and controller
    wind, sim, con, floris, floridyn, ta = setup(settings_file)
    # create settings struct
    set = Settings(wind, sim, con)
    wf, wind, sim, con, floris = prepareSimulation(set, wind, con, floridyn, floris, ta, sim)
    sim.n_sim_steps = 2
    wf = convert_wf_dict2windfarm(wf_dict)
    wf_old = deepcopy(wf)
    # M, wf = setUpTmpWFAndRun(set, wf, floris, wind)
end
nothing
