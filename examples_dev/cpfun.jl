# Copyright (c) 2025 Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

# Load DTU 10MW Cp and Ct tables (TSR vs blade pitch), provide interpolation and plot the results.
# Exposes:
#   cp(tsr, pitch_deg) -> Cp
#   ct(tsr, pitch_deg) -> Ct
#   cp_matrix, ct_matrix, tsr_values, pitch_values
using CSV, DataFrames
import DataFrames.SentinelArrays
using Interpolations
using ControlPlots, LaTeXStrings

# Load data
cp_df = CSV.read("data/DTU_10MW/cp.csv", DataFrame; header=1)
ct_df = CSV.read("data/DTU_10MW/ct.csv", DataFrame; header=1)

# First column name contains a slash; get it programmatically
col_ts = names(cp_df)[1]
_tsr_vals = Float64.(cp_df[!, col_ts])
_pitch_vals = parse.(Float64, string.(names(cp_df)[2:end]))

cp_matrix = Array{Float64}(undef, size(cp_df,1), size(cp_df,2)-1)
ct_matrix = Array{Float64}(undef, size(ct_df,1), size(ct_df,2)-1)
for (j,col) in enumerate(names(cp_df)[2:end])
	cp_matrix[:, j] = Float64.(cp_df[!, col])
	ct_matrix[:, j] = Float64.(ct_df[!, col])
end

# Create gridded linear interpolations
_cp_itp = interpolate((_tsr_vals, _pitch_vals), cp_matrix, Gridded(Linear()))
_ct_itp = interpolate((_tsr_vals, _pitch_vals), ct_matrix, Gridded(Linear()))

cp(tsr, pitch) = _cp_itp(tsr, pitch)
ct(tsr, pitch) = _ct_itp(tsr, pitch)

# Demo values
println("Cp(8.0, 0.75°) = ", cp(8.0, 0.75))
println("Ct(8.0, 0.75°) = ", ct(8.0, 0.75))

# Maxima (table resolution)
max_cp, lin_cp = findmax(cp_matrix)
idx_cp = CartesianIndices(cp_matrix)[lin_cp]
rcp, ccp = idx_cp.I
println("Max Cp = $(max_cp) at TSR=$( _tsr_vals[rcp] ), Pitch=$( _pitch_vals[ccp] )° (table grid)")
max_ct, lin_ct = findmax(ct_matrix)
idx_ct = CartesianIndices(ct_matrix)[lin_ct]
rct, cct = idx_ct.I
println("Max Ct = $(max_ct) at TSR=$( _tsr_vals[rct] ), Pitch=$( _pitch_vals[cct] )° (table grid)")

tsr_values::Vector{Float64}   = _tsr_vals
pitch_values::Vector{Float64} = _pitch_vals

# Build grids (no meshgrid in exported plt): result shape matches cp_matrix (length(tsr) × length(pitch))
TSR = repeat(tsr_values, 1, length(pitch_values))
PITCH = repeat(pitch_values', length(tsr_values), 1)
fig_cp = plt.figure("Cp Surface", figsize=(7,5))
ax_cp = fig_cp.add_subplot(1,1,1, projection="3d")
surf = ax_cp.plot_surface(TSR, PITCH, cp_matrix; cmap="viridis", linewidth=0, antialiased=true)
ax_cp.set_xlabel("TSR (-)")
ax_cp.set_ylabel("Pitch (deg)")
ax_cp.set_zlabel(L"C_p (-)")
ax_cp.set_title("DTU 10MW")
fig_cp.colorbar(surf, shrink=0.8, aspect=12, label=L"C_p", pad=0.1)  # pad increases distance from plot
plt.tight_layout()

# Ct surface plot (same style)
fig_ct = plt.figure("Ct Surface", figsize=(7,5))
ax_ct = fig_ct.add_subplot(1,1,1, projection="3d")
surf2 = ax_ct.plot_surface(TSR, PITCH, ct_matrix; cmap="viridis", linewidth=0, antialiased=true)
ax_ct.set_xlabel("TSR (-)")
ax_ct.set_ylabel("Pitch (deg)")
ax_ct.set_zlabel(L"C_t (-)")
ax_ct.set_title("DTU 10MW")
fig_ct.colorbar(surf2, shrink=0.8, aspect=12, label=L"C_t", pad=0.1)
plt.tight_layout()

settings_file, vis_file = get_default_project()[2:3]
# get the settings for the wind field, simulator, controller and turbines
_, _, _, _, _, _, tp = setup(settings_file)
@info "tp.cp_fun(8, 0.75) = $(tp.cp_fun(8, 0.75))"
