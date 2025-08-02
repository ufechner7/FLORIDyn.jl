# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using MAT
using Random
using Suppressor

if basename(pwd()) == "test"
    cd("..")
end

filename="test/randn.mat"
randn_vec=vec(matread(filename)["vec"])

mutable struct FileRNG <: AbstractRNG
    const data::Vector{Float64}
    idx::Int
end

FileRNG(data::Vector{Float64}) = FileRNG(data, 1)

function Base.rand(rng::FileRNG)
    x = rng.data[rng.idx]
    rng.idx += 1
    return x
end

function Base.rand(rng::FileRNG, ::Type{Float64})
    return rand(rng)
end

function Base.rand(rng::FileRNG, ::Type{wf}) where {wf}
    error("FileRNG only supports Float64")
end

Random.randn(rng::FileRNG) = rand(rng)

rng = FileRNG(randn_vec)
FLORIDyn.set_rng(rng)

# Define all available test files
all_test_files = [
    "test_dir.jl",
    "test_shear.jl", 
    "test_tit.jl",
    "test_vel.jl",
    "test_iterate.jl",
    "test_floridyn_cl.jl",
    "test_floris.jl",
    "test_correction.jl",
    "test_controller.jl",
    "test_settings.jl",
    "test_visualisation.jl",
    "test_parallel.jl",
    "test_copy_functions.jl",
    "aqua.jl"
]

# Files that need error suppression
suppress_error_files = [
    "test_dir.jl",
    "test_shear.jl",
    "test_tit.jl", 
    "test_vel.jl",
    "test_iterate.jl",
    "test_floridyn_cl.jl"
]

# Get test files to run from test_args
test_files_to_run = if isempty(ARGS)
    all_test_files
else
    # Filter to only include valid test files
    requested_files = String[]
    for arg in ARGS
        # Normalize the argument (add .jl if missing, ensure it's in test directory)
        normalized_arg = if endswith(arg, ".jl")
            arg
        else
            arg * ".jl"
        end
        
        # Check if it's a valid test file
        if normalized_arg in all_test_files
            push!(requested_files, normalized_arg)
        else
            @warn "Test file '$normalized_arg' not found. Available files: $(join(all_test_files, ", "))"
        end
    end
    
    if isempty(requested_files)
        @warn "No valid test files specified. Running all tests."
        all_test_files
    else
        requested_files
    end
end

@testset verbose=true "FLORIDyn Tests" begin
    # Run tests with error suppression for specific files
    files_needing_suppression = intersect(test_files_to_run, suppress_error_files)
    if !isempty(files_needing_suppression)
        @suppress_err begin
            for test_file in files_needing_suppression
                include(test_file)
            end
        end
    end
    
    # Run remaining tests without error suppression
    files_without_suppression = setdiff(test_files_to_run, suppress_error_files)
    for test_file in files_without_suppression
        include(test_file)
    end
end