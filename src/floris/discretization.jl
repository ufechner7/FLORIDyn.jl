# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    discretizeRotor(n_rp::Int) -> Tuple{Matrix{Float64}, Vector{Float64}}

Discretizes the rotor into a `n_rp` segments. The algorithm returns the normalized center location ∈ [-0.5, 0.5] and the
relative area the segment represents.

# Arguments
- `n_rp::Int`: The number of radial points to discretize the rotor into.

# Returns
- The tuple `(m_rp, w)` where:
  - `m_rp`: A matrix of size `(nC, 3)` where `nC` is the number of segments. The first column is all zeros, 
     the second and third columns contain the normalized radial positions.
  - `w`: A vector of weights corresponding to each segment, summing to approximately 1.

# Notes
- The algorithm returns the normalized center location in the range `[-0.5, 0.5]` and the
  relative area that each segment represents.
- The isocell algorithm is used, which may not yield exactly `n_rp` cells but aims to achieve a similar number.
- For details, see the publication by Masset et al.:
  [Masset et al. (2009)](https://orbi.uliege.be/bitstream/2268/91953/1/masset_isocell_orbi.pdf)
- The choice of `N1 = 3` is made here, but values of `4` or `5` are also viable options. The choice of `3` is close to optimal.
"""
function discretizeRotor(n_rp::Int)
    # DISCRETIZEROTOR discretizes the rotor plane into n_rp segments.
    #

    N1 = 3
    n = round(Int, sqrt(n_rp / N1))
    # eRP = n_rp - N1 * n^2  # Difference from calculated to actually applied

    # Radial thickness of each ring
    dltR = 1 / n
    nC = N1 * n^2

    # m_rp matrix: nC rows, 3 columns (initialize to zeros without creating a view)
    m_rp = Array{Float64}(undef, nC, 3)
    fill!(m_rp, 0.0)

    @inbounds for i in 1:n
        nR = (2 * i - 1) * N1
        # Sum of first i odd numbers = i^2, scaled by N1
        i_e = N1 * i^2
        i_s = i_e - nR + 1

        radius = 0.5 * dltR * (0.5 + (i - 1))
        step = 2π / nR
        for k in 0:(nR - 1)
            idx = i_s + k
            ang = (k + 1) * step
            m_rp[idx, 2] = radius * cos(ang)
            m_rp[idx, 3] = radius * sin(ang)
        end
  end

  w = Vector{Float64}(undef, nC)
  fill!(w, 1 / nC)

    return m_rp, w
end
