# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    discretizeRotor(nRP::Int)

Discretizes the rotor into a `nRP` segments. The algorithm returns the normalized center location ∈ [-0.5, 0.5] and the
relative area the segment represents.

# Arguments
- `nRP::Int`: The number of radial points to discretize the rotor into.

# Returns
- The tuple `(RPs, w)` where:
  - `RPs`: A matrix of size `(nC, 3)` where `nC` is the number of segments. The first column is all zeros, 
     the second and third columns contain the normalized radial positions.
  - `w`: A vector of weights corresponding to each segment, summing to approximately 1.

# Notes
- The algorithm returns the normalized center location in the range `[-0.5, 0.5]` and the
  relative area that each segment represents.
- The isocell algorithm is used, which may not yield exactly `nRP` cells but aims to achieve a similar number.
- For details, see the publication by Masset et al.:
  [Masset et al. (2009)](https://orbi.uliege.be/bitstream/2268/91953/1/masset_isocell_orbi.pdf)
- The choice of `N1 = 3` is made here, but values of `4` or `5` are also viable options. The choice of `3` is close to optimal.
"""
function discretizeRotor(nRP::Int)
    # DISCRETIZEROTOR discretizes the rotor plane into nRP segments.
    #

    N1 = 3
    n = round(Int, sqrt(nRP / N1))
    # eRP = nRP - N1 * n^2  # Difference from calculated to actually applied

    # Radial thickness of each ring
    dltR = 1 / n
    nC = N1 * n^2

    # RPs matrix: nC rows, 3 columns
    # Columns:
    # RPs[:, 1] corresponds to RPs(:,1) in MATLAB, but unused in original code (remains zeros)
    RPs = zeros(nC, 3)

    for i in 1:n
        nR = (2 * i - 1) * N1
        i_e = sum(((2 .* (1:i) .- 1) .* N1))
        i_s = i_e - nR + 1

        phi = (1:nR) .* (2 * π) ./ nR
        RPs[i_s:i_e, 2] .= 0.5 .* cos.(phi) .* dltR .* (0.5 + (i - 1))
        RPs[i_s:i_e, 3] .= 0.5 .* sin.(phi) .* dltR .* (0.5 + (i - 1))
    end

    w = fill(1 / nC, nC)

    return RPs, w
end
