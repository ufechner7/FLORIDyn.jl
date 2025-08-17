# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    IterateOPs_model

Abstract supertype for all operational point iteration algorithms in FLORIDyn.

This abstract type defines the interface for different strategies used to advance 
operational points through the wind field during time-stepping simulations. All 
concrete iteration models must be subtypes of this abstract type.

# Purpose
The iteration models determine how operational points (OPs) move through space 
and time, affecting:
- Wake propagation dynamics
- Spatial discretization accuracy
- Computational efficiency
- Physical representation of wind farm interactions

# Implementation
Concrete subtypes implement specific iteration strategies through method dispatch 
on functions like [`iterateOPs!`](@ref). Each model represents a different 
approach to handling the temporal and spatial evolution of operational points.

# Available Models
- [`IterateOPs_basic`](@ref): Basic time-stepping with simple advection
- [`IterateOPs_average`](@ref): Averaged dynamics for stability
- [`IterateOPs_buffer`](@ref): Buffered approach for memory efficiency
- [`IterateOPs_maximum`](@ref): Maximum value-based iteration
- [`IterateOPs_weighted`](@ref): Weighted interpolation method
"""
abstract type IterateOPs_model end

"""
    IterateOPs_average <: IterateOPs_model

Operational point iteration model using averaged dynamics.

This iteration strategy employs averaging techniques to advance operational 
points through the wind field, providing enhanced numerical stability and 
smoother wake evolution compared to basic methods.

# Algorithm Characteristics
- **Stability**: Improved numerical stability through temporal averaging
- **Smoothness**: Reduces oscillations in wake dynamics
- **Computational Cost**: Moderate overhead due to averaging operations
- **Accuracy**: Good balance between stability and physical representation

# Use Cases
Recommended for:
- Simulations requiring smooth wake evolution
- Cases with high turbulence or complex wind conditions
- Long-duration simulations where stability is critical
- Research applications focusing on ensemble statistics

# Mathematical Approach
The averaging process involves temporal or spatial averaging of 
relevant quantities (velocities, deflections, turbulence) before applying 
the advancement step, resulting in more stable operational point trajectories.

# Notes
This model may require additional computational resources compared to basic 
methods but provides better stability characteristics for challenging 
simulation scenarios.
"""
struct IterateOPs_average <: IterateOPs_model end

"""
    IterateOPs_basic <: IterateOPs_model

Basic operational point iteration model with simple time-stepping.

This is the fundamental iteration strategy that advances operational points 
using direct time-stepping based on local wind velocities and wake deflection 
effects. It provides the core functionality for FLORIDyn simulations.

# Algorithm Characteristics
- **Simplicity**: Straightforward implementation with minimal overhead
- **Performance**: Fastest execution among available iteration models
- **Accuracy**: Direct physical representation of wake advection
- **Memory**: Minimal memory requirements

# Implementation
The basic algorithm performs:
1. Downwind advection based on local wind speed
2. Crosswind deflection using wake centerline calculations
3. Coordinate transformation to world coordinates
4. Temporal advancement through circular shifting
5. Spatial reordering to maintain downstream position order

# Use Cases
Recommended for:
- Standard wind farm simulations
- Performance-critical applications
- Validation studies against reference data
- Initial model development and testing

# Mathematical Foundation
Uses explicit time-stepping with:
```
Δx = U × Δt × advection_factor
```
where operational points move downstream based on local wind conditions.

# Notes
This model serves as the reference implementation and baseline for 
comparison with other iteration strategies. See [`iterateOPs!`](@ref)
for detailed implementation.
"""
struct IterateOPs_basic <: IterateOPs_model end

"""
    IterateOPs_buffer <: IterateOPs_model

Operational point iteration model with buffered memory management.

This iteration strategy implements buffering techniques to optimize memory 
usage and computational efficiency during operational point advancement, 
particularly beneficial for large-scale wind farm simulations.

# Algorithm Characteristics
- **Memory Efficiency**: Optimized memory access patterns
- **Scalability**: Better performance for large numbers of operational points
- **Caching**: Intelligent buffering of frequently accessed data
- **Computational Cost**: Reduced overhead for memory-intensive operations

# Use Cases
Recommended for:
- Large wind farms with many turbines
- Memory-constrained computing environments
- High-resolution simulations with dense operational point grids
- Production simulations requiring optimal resource utilization

# Implementation Strategy
The buffering approach manages operational point data through:
- Efficient memory allocation patterns
- Reduced data copying operations
- Optimized access to state matrices
- Strategic caching of intermediate results

# Performance Benefits
- Improved cache locality for better CPU performance
- Reduced memory bandwidth requirements
- Better scaling with problem size
- Lower memory fragmentation

# Notes
This model is particularly effective when computational resources are 
limited or when dealing with very large simulation domains.
"""
struct IterateOPs_buffer <: IterateOPs_model end

"""
    IterateOPs_maximum <: IterateOPs_model

Operational point iteration model using maximum value selection.

This iteration strategy employs maximum value-based decision making during 
operational point advancement, potentially useful for conservative estimates 
or worst-case scenario analysis in wind farm simulations.

# Algorithm Characteristics
- **Conservative Approach**: Tends toward maximum/conservative values
- **Robustness**: Provides bounds on simulation behavior
- **Special Cases**: Handles extreme conditions effectively
- **Analysis**: Useful for sensitivity and worst-case studies

# Use Cases
Recommended for:
- Conservative design analysis
- Worst-case scenario evaluation
- Risk assessment studies
- Validation of simulation bounds
- Research on extreme wind conditions

# Mathematical Approach
The maximum selection process applies maximum operators to relevant 
quantities during the iteration step, which may include:
- Maximum wind speeds in the vicinity
- Maximum deflection values
- Maximum turbulence intensities
- Conservative time step selection

# Applications
Particularly useful in:
- Safety factor determination
- Conservative power estimation
- Extreme load analysis
- Uncertainty quantification studies

# Notes
This model may produce more conservative results compared to other 
iteration strategies and should be used when understanding bounds 
on simulation behavior is important.
"""
struct IterateOPs_maximum <: IterateOPs_model end

"""
    IterateOPs_weighted <: IterateOPs_model

Operational point iteration model using weighted interpolation.

This iteration strategy employs sophisticated weighted interpolation 
techniques to advance operational points, providing enhanced accuracy 
through spatial and temporal weighting of relevant physical quantities.

# Algorithm Characteristics
- **High Accuracy**: Superior interpolation accuracy
- **Smoothness**: Smooth transitions between operational points
- **Computational Cost**: Higher due to interpolation calculations
- **Flexibility**: Adaptable weighting schemes

# Use Cases
Recommended for:
- High-accuracy research simulations
- Detailed wake interaction studies
- Validation against experimental data
- Applications requiring smooth field representations
- Advanced control algorithm development

# Interpolation Strategy
The weighted approach typically involves:
- Distance-based spatial weighting
- Temporal interpolation of quantities
- Multi-point interpolation schemes
- Adaptive weight calculation based on local conditions

# Mathematical Foundation
Uses interpolation weights w_i such that:
```
quantity_interpolated = Σ w_i × quantity_i
```
where weights satisfy Σ w_i = 1 and are computed based on distance, 
time, or other relevant metrics.

# Accuracy Benefits
- Reduced numerical diffusion
- Better preservation of wake structures
- Improved representation of turbulence effects
- Enhanced spatial resolution

# Notes
This model provides the highest accuracy among available iteration 
strategies but requires additional computational resources. Best suited 
for applications where accuracy is prioritized over computational speed.
"""
struct IterateOPs_weighted <: IterateOPs_model end

"""
    FLORISBuffers

Pre-allocated buffers for the runFLORIS computation to minimize allocations.

This struct also persists result arrays so callers can read outputs without
allocations. After calling `runFLORIS`, the following fields contain results:

- `T_red_arr::Vector{Float64}`: Per-turbine velocity reduction factors. For a
    single-turbine run, length is 1 and `T_red_arr[1]` is the scalar reduction.
- `T_aTI_arr::Vector{Float64}`: Added turbulence intensity from upstream wakes.
    For N turbines, length is `max(N-1, 0)`. Empty for single-turbine runs.
- `T_Ueff::Vector{Float64}`: Effective wind speed at the last turbine as a
    length-1 vector (multi-turbine case). Empty for single-turbine runs.
- `T_weight::Vector{Float64}`: Gaussian weight factors used for wake overlap.
    For N turbines, length is `max(N-1, 0)`. Empty for single-turbine runs.
"""
mutable struct FLORISBuffers
    tmp_RPs::Matrix{Float64}
    rotor_pts::Matrix{Float64}
    # Preallocated arrays for getVars! outputs
    sig_y::Vector{Float64}
    sig_z::Vector{Float64}
    x_0::Vector{Float64}
    delta::Matrix{Float64}   # n×2
    pc_y::Vector{Float64}
    pc_z::Vector{Float64}
    cw_y::Vector{Float64}
    cw_z::Vector{Float64}
    phi_cw::Vector{Float64}
    r_cw::Vector{Float64}
    core::Vector{Bool}
    nw::Vector{Bool}
    fw::Vector{Bool}
    tmp_RPs_r::Vector{Float64}
    gaussAbs::Vector{Float64}
    gaussWght::Vector{Float64}
    exp_y::Vector{Float64}
    exp_z::Vector{Float64}
    not_core::Vector{Bool}
    # Result arrays (persisted in buffers to avoid fresh allocations)
    T_red_arr::Vector{Float64}
    T_aTI_arr::Vector{Float64}
    T_Ueff::Vector{Float64}    # length 1 when set
    T_weight::Vector{Float64}
end


"""
    UnifiedBuffers

Unified buffer struct containing all arrays needed by interpolateOPs! and setUpTmpWFAndRun!.

# Fields
- `dist_buffer::Vector{Float64}`: Distance calculations for interpolateOPs!
- `sorted_indices_buffer::Vector{Int}`: Sorted indices for interpolateOPs!
- `M_buffer::Matrix{Float64}`: Main result buffer for setUpTmpWFAndRun!
- `iTWFState_buffer::Vector{Float64}`: Turbine wind field state buffer
- `tmp_Tpos_buffer::Matrix{Float64}`: Temporary turbine position buffer
- `tmp_WF_buffer::Matrix{Float64}`: Temporary wind field buffer
- `tmp_Tst_buffer::Matrix{Float64}`: Temporary turbine state buffer
- `dists_buffer::Vector{Float64}`: Distance buffer for setUpTmpWFAndRun!
- `plot_WF_buffer::Matrix{Float64}`: Wind field plotting buffer
- `plot_OP_buffer::Matrix{Float64}`: Observation point plotting buffer
- `floris_buffers`: Pre-allocated FLORIS computation buffers for wind wake calculations
"""
struct UnifiedBuffers
    dist_buffer::Vector{Float64}
    sorted_indices_buffer::Vector{Int}
    M_buffer::Matrix{Float64}
    iTWFState_buffer::Vector{Float64}
    tmp_Tpos_buffer::Matrix{Float64}
    tmp_WF_buffer::Matrix{Float64}
    tmp_Tst_buffer::Matrix{Float64}
    dists_buffer::Vector{Float64}
    plot_WF_buffer::Matrix{Float64}
    plot_OP_buffer::Matrix{Float64}
    floris_buffers::FLORISBuffers  # Will be FLORISBuffers when FLORIS is available
    gp::Any              # Optional WindFarm buffer for grid-point computations
end

