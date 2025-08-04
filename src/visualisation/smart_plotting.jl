# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

if Threads.nthreads() > 1
    include("remote_plotting.jl") 
    init_plotting()  # This sets up workers and remote plotting capabilities   
end

