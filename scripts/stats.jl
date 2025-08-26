# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# List all methods of FLORIDyn
using FLORIDyn

exported_names = names(FLORIDyn)
exported_functions = filter(n -> isdefined(FLORIDyn, n) && isa(getfield(FLORIDyn, n), Function), exported_names)
println("Exported methods:\n")
total = 0
for fun in exported_functions
    global total
    f = getfield(FLORIDyn, fun)
    mes = methods(f)
    for me in mes
        println(replace(split(repr(me),'@')[1], "FLORIDyn." => ""))
        total += 1
    end
end
println("\nTotal: $total")

function check_exported_docs(mod::Module)
    exported_symbols = names(mod, all=false)
    doc_status = Dict{Symbol,Bool}()
    for sym in exported_symbols
        #if isa(getfield(mod, sym), Function)
            doc_status[sym] = Base.Docs.hasdoc(mod, sym)
        #end
    end
    return doc_status
end

# Usage example:
results = check_exported_docs(FLORIDyn)
undocumented = filter(kv -> !kv[2], results)
println("\nUndocumented exported symbols: \n", keys(undocumented))
