# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# build and display the html documentation locally
# you must have installed the package LiveServer in your global environment
using Pkg

function globaldependencies()
    projectpath = Pkg.project().path
    basepath, _ = splitdir(projectpath)
    Pkg.activate()
    globaldependencies = keys(Pkg.project().dependencies)
    Pkg.activate(basepath)
    globaldependencies
end

if !("LiveServer" in globaldependencies())
    println("Installing LiveServer globally!")
    run(`julia -e 'using Pkg; Pkg.add("LiveServer")'`)
end

if !("Documenter" ∈ keys(Pkg.project().dependencies))
    using TestEnv
    TestEnv.activate()
end
using LiveServer; servedocs(launch_browser=true)
