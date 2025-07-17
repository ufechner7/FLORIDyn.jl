# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

function discretizeRotor(nRP)
    # DISCRETIZEROTOR discretizes the rotor plane into nRP segments.
    # The algorithm returns the normalized center location ∈ [-0.5, 0.5] and the
    # relative area the segment represents.
    #
    # The isocell algorithm
    # It faces certain limitations and aims to achieve nRP cells but might end
    # up with a slightly different number. For details, see the publication of
    # Masset et al.:
    #  https://orbi.uliege.be/bitstream/2268/91953/1/masset_isocell_orbi.pdf
    # We choose N1 = 3 here, 4 or 5 are also viable options, 3 is close to
    # optimal

    N1 = 3
    n = round(Int, sqrt(nRP / N1))
    # eRP = nRP - N1 * n^2  # Difference from calculated to actually applied

    # Radial thickness of each ring
    dltR = 1 / n
    nC = N1 * n^2

    # RPs matrix: nC rows, 3 columns
    # Columns: (I think MATLAB first is 1-based; in Julia too)
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
