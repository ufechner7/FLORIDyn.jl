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
    authors="Uwe Fechner <uwe.fechner.msc@gmail.com>, Marcus Becker and contributors",
    repo="https://github.com/ufechner7/FLORIDyn.jl/blob/{commit}{path}#{line}",
    sitename="FLORIDyn.jl",
    checkdocs=:none,
    format=Documenter.HTML(;
      assets = String["assets/custom.css"],
      canonical = "https://ufechner7.github.io/FLORIDyn.jl/dev/",
      repolink = "https://github.com/ufechner7/FLORIDyn.jl",
      prettyurls=get(ENV, "CI", "false") == "true",
    ),
    pages=[
        "Home" => "index.md",
        "Exported Types" => "types.md",
        "Low-Level Functions" => "functions.md",
        "High-Level Functions" => "functions_II.md",
        "Settings" => "settings.md",
        "Video Creation" => "video_creation.md", 
        "User Guide" => "user_guide.md",
    "Examples" => "examples.md",
        "Developer notes" => "developer.md",
        "FAQ" => "faq.md"
    ],
)

deploydocs(;
    repo="github.com/ufechner7/FLORIDyn.jl",
    devbranch="main",
    push_preview=true,
)
