# Load DTU 10MW Cp and Ct tables (TSR vs blade pitch) and provide interpolation.
# Exposes:
#   cp(tsr, pitch_deg) -> Cp
#   ct(tsr, pitch_deg) -> Ct
#   cp_matrix, ct_matrix, tsr_values, pitch_values
using CSV, DataFrames
using Interpolations

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

const tsr_values = _tsr_vals
const pitch_values = _pitch_vals
