# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using MAT
using Random

if basename(pwd()) == "test"
    cd("..")
end

filename="test/randn.mat"
randn_vec=vec(matread(filename)["vec"])

# const DNG = Iterators.Stateful(randn_vec)
# randn() = popfirst!(DNG)
# randn(n) = collect(Iterators.take(DNG, n))

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

function Base.rand(rng::FileRNG, ::Type{T}) where {T}
    error("FileRNG only supports Float64")
end

# Optionally, support randn
Random.randn(rng::FileRNG) = rand(rng)

rng = FileRNG(randn_vec)


# include("test_windfield.jl")