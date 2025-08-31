# Load DTU 10MW Cp and Ct tables (TSR vs blade pitch) and provide interpolation.
# Exposes:
#   cp(tsr, pitch_deg) -> Cp
#   ct(tsr, pitch_deg) -> Ct
#   cp_matrix, ct_matrix, tsr_values, pitch_values
using CSV, DataFrames
import DataFrames.SentinelArrays
using Interpolations
using ControlPlots

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
const _cp_itp = interpolate((_tsr_vals, _pitch_vals), cp_matrix, Gridded(Linear()))
const _ct_itp = interpolate((_tsr_vals, _pitch_vals), ct_matrix, Gridded(Linear()))

cp(tsr, pitch) = _cp_itp(tsr, pitch)
ct(tsr, pitch) = _ct_itp(tsr, pitch)

# Demo values
println("Cp(5.0, 10.0°) = ", cp(5.0, 10.0))
println("Ct(5.0, 10.0°) = ", ct(5.0, 10.0))

tsr_values::Vector{Float64}   = _tsr_vals
pitch_values::Vector{Float64} = _pitch_vals

# Build grids (no meshgrid in exported plt): result shape matches cp_matrix (length(tsr) × length(pitch))
TSR = repeat(tsr_values, 1, length(pitch_values))
PITCH = repeat(pitch_values', length(tsr_values), 1)
fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(1,1,1, projection="3d")
surf = ax.plot_surface(TSR, PITCH, cp_matrix; cmap="viridis", linewidth=0, antialiased=true)
ax.set_xlabel("TSR (-)")
ax.set_ylabel("Pitch (deg)")
ax.set_zlabel("Cp (-)")
ax.set_title("DTU 10MW Cp Surface")
fig.colorbar(surf, shrink=0.8, aspect=12, label="Cp", pad=0.1)  # pad increases distance from plot
plt.tight_layout()
