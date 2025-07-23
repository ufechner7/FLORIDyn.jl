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

@testset verbose=true "FLORIDyn Tests" begin
    @suppress_err begin
        include("test_dir.jl")
        include("test_shear.jl")
        include("test_tit.jl")
        include("test_vel.jl")
    end
    include("test_floris.jl")
    include("test_correction.jl")
    include("test_floridyn_cl.jl")
    include("test_init.jl")
end