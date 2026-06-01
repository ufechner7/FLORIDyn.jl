# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Example of a simple Julia script for local data processing.

identity_func(x) = x

function process_local(data)
    identity_func(data)
end

data = rand(UInt8, 1024^2)  # 1MB random data

# Pre-compile by calling once
process_local(data)

# Measure local processing time
@time process_local(data)
nothing
