# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Determine if multiple Y series are provided and check lengths for three different parameter styles.

function dummy_plot(X, Ys...; xlabel="", ylabel="", ylabels=nothing, labels=nothing, xlims=nothing, ylims=nothing, ann=nothing, 
    scatter=false, title="", fig="", ysize=14, pltctrl=nothing)
    if length(Ys) > 1
        println("Multiple Y series detected: ", length(Ys))
        for (i, Y) in pairs(Ys)
            println("Length of X ($(length(X))),  length of Ys[$i] ($(length(Y))).")
            if length(X) != length(Y)
                throw(ArgumentError("Length of X ($(length(X))) does not match length of Ys[$i] ($(length(Y)))."))
            end
        end
    else
        println("Single Y series detected.")
        if length(Ys) == 1 && !isa(Ys[1], AbstractArray{<:AbstractArray})
            Y = Ys[1]
            if length(X) != length(Y)
                throw(ArgumentError("Length of X ($(length(X))) does not match length of Ys[1] ($(length(Y)))."))
            end
        else
            for (i, Y) in pairs(Ys[1])
                println("Length of X ($(length(X))),  length of Ys[1][$i] ($(length(Y))).")
                if length(X) != length(Y)
                    throw(ArgumentError("Length of X ($(length(X))) does not match length of Ys[1][$i] ($(length(Y)))."))
                end
            end
        end
    end

end

# Create test data
X = [1.0, 2.0, 3.0, 4.0, 5.0]
Y1 = [2.0, 4.0, 3.0, 5.0, 4.5]
Y2 = [1.0, 3.0, 2.0, 4.0, 3.5]

# Test basic functionality with single Y series (should work)
dummy_plot(X, Y1; xlabel="X Values", title="Test Plot", fig="Basic Test")

# Test functionality with multiple Y series (should work)
dummy_plot(X, Y1, Y2; xlabel="X Values", title="Test Plot", fig="Basic Test")

dummy_plot(X, [Y1, Y2]; xlabel="X Values", title="Test Plot", fig="Basic Test")
