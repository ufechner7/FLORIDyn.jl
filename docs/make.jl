# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

using FLORIDyn
using Pkg
if ("TestEnv" ∈ keys(Pkg.project().dependencies))
    if ! ("Documenter" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
end
using Documenter

DocMeta.setdocmeta!(FLORIDyn, :DocTestSetup, :(using FLORIDyn); recursive=true)

makedocs(;
    modules=[FLORIDyn],
    authors="Uwe Fechner <uwe.fechner.msc@gmail.com>, Markus Becker and contributors",
    repo="https://github.com/ufechner7/FLORIDyn.jl/blob/{commit}{path}#{line}",
    sitename="FLORIDyn.jl",
    checkdocs=:none,
    format=Documenter.HTML(;
      repolink = "https://github.com/ufechner7/FLORIDyn.jl",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ufechner7.github.io/FLORIDyn.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Exported Types" => "types.md",
        "Exported Functions" => "functions.md",
    ],
)

deploydocs(;
    repo="github.com/ufechner7/FLORIDyn.jl",
    devbranch="main",
)
