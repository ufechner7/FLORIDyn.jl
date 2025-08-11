# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using Aqua

@testset verbose=true "AQUA Quality Assurance Tests" begin
    Aqua.test_all(FLORIDyn)
end