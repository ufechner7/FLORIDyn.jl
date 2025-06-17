# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Test
using LinearAlgebra
using Random
using MAT

Random.seed!(1234)
filename="test/randn.mat"
randn_vec=vec(matread(filename)["vec"])

include("test_windfield.jl")