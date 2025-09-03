# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# This script shows how to plot multiple time series in one plot

using ControlPlots, LaTeXStrings

X = 0:0.1:2pi
Y1 = sin.(X)
Y2 = cos.(X)
p = plot(X, [Y1, Y2], xlabel=L"\alpha = [0..2\pi]", ylabel="Function Value", labels=["sin","cos"], fig="Dual")