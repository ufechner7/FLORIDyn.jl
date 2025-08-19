#!/usr/bin/env julia
# Simple launcher for the FLORIDyn project selector
import Pkg
try
    # Prefer running in project context if available
    Base.active_project() === nothing && Pkg.activate(@__DIR__)
catch
end
using FLORIDyn
select_project()
