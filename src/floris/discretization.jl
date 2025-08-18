# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    discretizeRotor(n_rp::Int) -> Tuple{Matrix{Float64}, Vector{Float64}}

Discretizes the rotor into `n_rp` segments using the isocell algorithm.

Memoization: results are cached per thread. Repeated calls with the same `n_rp`
on the same thread reuse the cached arrays (no lock needed). Do not mutate the
returned arrays, as they are shared within the thread.

# Arguments
- `n_rp::Int`: The number of radial points to discretize the rotor into.

# Returns
- `(m_rp, w)` where:
  - `m_rp::Matrix{Float64}`: Size `(nC, 3)`; first column zeros, columns 2–3 are
    normalized coordinates in `[-0.5, 0.5]`.
  - `w::Vector{Float64}`: Weights per cell that sum to approximately 1.

# Notes
- Per-thread cache avoids contention; different threads may compute and hold
  their own cached copies for the same `n_rp`.
- The isocell algorithm may not yield exactly `n_rp` cells but aims for a similar number.
- For details, see: Masset et al. (2009)
  https://orbi.uliege.be/bitstream/2268/91953/1/masset_isocell_orbi.pdf
- The choice `N1 = 3` is used here; values of 4 or 5 are also viable.
"""
function discretizeRotor(n_rp::Int)
  cache = _get_discretizeRotor_thread_cache()
  return get!(cache, n_rp) do
    _compute_discretizeRotor(n_rp)
  end
end

# Per-thread memoization caches using a growable vector of Dicts stored in a Ref.
# We guard only the growth operation with a lock; steady-state lookups are lock-free.
const _discretizeRotor_caches_ref = Base.RefValue{Vector{Dict{Int, Tuple{Matrix{Float64}, Vector{Float64}}}}}(Vector{Dict{Int, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, 0))
const _discretizeRotor_cache_lock = ReentrantLock()

@inline function _get_discretizeRotor_thread_cache()
  caches = _discretizeRotor_caches_ref[]
  tid = Base.Threads.threadid()
  if length(caches) < tid
    lock(_discretizeRotor_cache_lock) do
      caches2 = _discretizeRotor_caches_ref[]  # re-read in case another thread updated
      if length(caches2) < tid
        # Grow to at least current nthreads(), preserving previous entries
        newlen = max(Base.Threads.nthreads(), tid)
        newc = Vector{Dict{Int, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, newlen)
        @inbounds for i in 1:newlen
          newc[i] = i <= length(caches2) ? caches2[i] : Dict{Int, Tuple{Matrix{Float64}, Vector{Float64}}}()
        end
        _discretizeRotor_caches_ref[] = newc
        caches = newc
      else
        caches = caches2
      end
    end
  end
  return caches[tid]
end

@inline function _compute_discretizeRotor(n_rp::Int)
  # DISCRETIZEROTOR discretizes the rotor plane into n_rp segments.

  N1 = 3
  n = round(Int, sqrt(n_rp / N1))

  # Radial thickness of each ring
  dltR = 1 / n
  nC = N1 * n^2

  # m_rp matrix: nC rows, 3 columns (initialize to zeros)
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
