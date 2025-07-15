# list all functions of FLORIDyn
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
        println(split(repr(me),'@')[1])
        total += 1
    end
end
println("\nTotal: $total")