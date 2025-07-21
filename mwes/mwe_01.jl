# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

A =  [1500.0, 2400.0]
B =  ones(200, 2) * 1000.0
C = A' .- B .+ 100*rand(200)

 distOP_iiT = sum((C).^2, dims=2)
 I_op = argmin(distOP_iiT)