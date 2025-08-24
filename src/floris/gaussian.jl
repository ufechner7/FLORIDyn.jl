# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

#=
This file contains core functions for the Gaussian wake model implementation in FLORIS.

Functions and structures defined in this file:
- calcCt: Calculate thrust coefficient from axial induction factor (scalar and vectorized versions)
- FLORISBuffers: Constructor for pre-allocated computation buffers
- getVars!: Compute Gaussian wake widths, deflection, potential-core radii, and onset distance
- centerline!: Compute cross-wind wake deflection at observation points
- States (struct): Container for state variable names and counts 
- States(): Default constructor for States struct with predefined variables
- init_states: Initialize state arrays for wind farm simulation
- getPower: Calculate power output of wind turbines accounting for yaw effects
- getUadv: Calculate advection speed factor based on Zong & Porté-Agel 2020

This file provides the mathematical foundation for the Gaussian wake model, including
wake expansion calculations, deflection modeling, and state management. The main
runFLORIS! function and its helpers have been moved to runfloris.jl for better
code organization.
=#

"""
    calcCt(a::Number, _) -> Float64
    calcCt(a::AbstractVector, _) -> AbstractVector

Calculate the thrust coefficient ct = 4a(1-a).
"""
@inline function calcCt(a::Number, _)::Float64
    return 4 * a * (1 - a)
end

@inline function calcCt(a::AbstractVector, _)
    return 4 .* a .* (1 .- a)
end


function FLORISBuffers(n_pts::Int)
    return FLORISBuffers(
        Matrix{Float64}(undef, n_pts, 3),  # tmp_RPs
        Matrix{Float64}(undef, n_pts, 3),  # rotor_pts
        Vector{Float64}(undef, n_pts),     # sig_y
        Vector{Float64}(undef, n_pts),     # sig_z
        Vector{Float64}(undef, n_pts),     # x_0
        Matrix{Float64}(undef, n_pts, 2),  # delta
        Vector{Float64}(undef, n_pts),     # pc_y
        Vector{Float64}(undef, n_pts),     # pc_z
        Vector{Float64}(undef, n_pts),     # cw_y
        Vector{Float64}(undef, n_pts),     # cw_z
        Vector{Float64}(undef, n_pts),     # phi_cw
        Vector{Float64}(undef, n_pts),     # r_cw
        Vector{Bool}(undef, n_pts),        # core
        Vector{Bool}(undef, n_pts),        # nw
        Vector{Bool}(undef, n_pts),        # fw
        Vector{Float64}(undef, n_pts),     # tmp_RPs_r
        Vector{Float64}(undef, n_pts),     # gaussAbs
        Vector{Float64}(undef, n_pts),     # gaussWght
        Vector{Float64}(undef, n_pts),     # exp_y
        Vector{Float64}(undef, n_pts),     # exp_z
        Vector{Bool}(undef, n_pts),        # not_core
        Float64[],                         # T_red_arr (size set per call)
        Float64[],                         # T_aTI_arr (size set per call)
        Float64[],                         # T_Ueff (size 0 or 1)
        Float64[],                         # T_weight (size set per call)
    )
end

"""
        getVars!(sig_y, sig_z, x_0, delta, pc_y, pc_z, rps, c_t, yaw, ti, ti0, floris::Floris, d_rotor)

Compute Gaussian wake widths, deflection, potential-core radii, and onset distance at observation points, in-place.

# Output Arguments
- `sig_y::AbstractVector{<:Real}` (length n): Lateral Gaussian width σ_y at each point [m]
- `sig_z::AbstractVector{<:Real}` (length n): Vertical Gaussian width σ_z at each point [m]
- `x_0::AbstractVector{<:Real}` (length n): Onset distance of the far-wake x₀ [m]
- `delta::AbstractMatrix{<:Real}` (length n×2): Deflection components `[Δy, Δz]` [m]
- `pc_y::AbstractVector{<:Real}` (length n): Potential-core radius in y at each point [m]
- `pc_z::AbstractVector{<:Real}` (length n): Potential-core radius in z at each point [m]

# Input Arguments
- `rps::AbstractMatrix` (n×3): Observation points in wake-aligned frame; columns are `[x_downstream, y_cross, z_cross]` [m]
- `c_t::Union{Number,AbstractVector}`: Thrust coefficient Ct (scalar or length n) [-]
- `yaw::Union{Number,AbstractVector}`: Yaw misalignment (scalar or length n) [rad]
- `ti::Union{Number,AbstractVector}`: Local turbulence intensity TI at turbine (scalar or length n) [-]
- `ti0::Union{Number,AbstractVector}`: Ambient turbulence intensity TI₀ (scalar or length n) [-]
- `floris::Floris`: FLORIS Gaussian model parameters; see [`Floris`](@ref)
- `d_rotor::Real`: Rotor diameter D [m]

# Behavior
- Supports scalar or per-point values for `c_t`, `yaw`, `ti`, `ti0`; scalars are broadcast to all points.
- Uses `floris.k_a`, `floris.k_b`, `floris.alpha`, `floris.beta` to compute per-point
    `x₀`, `σ_y`, `σ_z`, deflection `Δy` (here `Δz` is set to 0), and potential-core radii `pc_y`, `pc_z`.
- No heap allocations beyond the provided outputs; results are written in-place and the function returns `nothing`.

# Notes
- Only `rps[:, 1]` (downstream distance) is used by this implementation; `rps[:, 2:3]` are ignored.
- `delta` must have at least 2 columns; only columns 1:2 are written.
- Units: distances in meters, angles in radians, intensities and `Ct` are dimensionless.

# Example
```julia
n = size(RPs, 1)
sig_y = similar(RPs[:, 1])
sig_z = similar(RPs[:, 1])
x0    = similar(RPs[:, 1])
delta = zeros(n, 2)
pc_y  = similar(RPs[:, 1])
pc_z  = similar(RPs[:, 1])
getVars!(sig_y, sig_z, x0, delta, pc_y, pc_z, RPs, Ct, yaw, TI, TI0, floris, D)
```

Returns
- `nothing` — all results are written into the provided arrays.
"""
function getVars!(sig_y::AbstractVector{<:Real},
                  sig_z::AbstractVector{<:Real},
                  x_0::AbstractVector{<:Real},
                  delta::AbstractMatrix{<:Real},
                  pc_y::AbstractVector{<:Real},
                  pc_z::AbstractVector{<:Real},
                  rps::Union{Matrix, Adjoint, SubArray}, c_t, yaw, ti, ti0,
                  floris::Floris, d_rotor)
    # Parameters and constants
    k_a   = floris.k_a
    k_b   = floris.k_b
    alpha = floris.alpha
    beta  = floris.beta
    invsqrt8 = 1 / sqrt(8)
    sqrt2 = sqrt(2)

    n = size(rps, 1)
    @assert length(sig_y) >= n && length(sig_z) >= n && length(x_0) >= n &&
            size(delta, 1) >= n && size(delta, 2) >= 2 && length(pc_y) >= n && length(pc_z) >= n

    is_scalar_ct  = isa(c_t, Number)
    is_scalar_yaw = isa(yaw, Number)
    is_scalar_ti  = isa(ti, Number)
    is_scalar_ti0 = isa(ti0, Number)

    # Prepare x_0: scalar fill when inputs are scalars, else per-point
    if is_scalar_ct && is_scalar_yaw && is_scalar_ti && is_scalar_ti0
        Ct = c_t; yw = yaw
        I = sqrt(ti^2 + ti0^2)
        cosyaw = cos(yw)
        denom = sqrt2 * (alpha * I + beta * (1 - sqrt(1 - Ct)))
        x0 = (cosyaw * (1 + sqrt(1 - Ct))) / denom * d_rotor
        @inbounds fill!(x_0, x0)
    else
        @inbounds for i in 1:n
            Ct_i  = is_scalar_ct  ? c_t  : c_t[i]
            yaw_i = is_scalar_yaw ? yaw : yaw[i]
            ti_i  = is_scalar_ti  ? ti  : ti[i]
            ti0_i = is_scalar_ti0 ? ti0 : ti0[i]
            I_i = sqrt(ti_i^2 + ti0_i^2)
            cosyaw_i = cos(yaw_i)
            denom = sqrt2 * (alpha * I_i + beta * (1 - sqrt(1 - Ct_i)))
            x_0[i] = (cosyaw_i * (1 + sqrt(1 - Ct_i))) / denom * d_rotor
        end
    end

    # Main loop
    @inbounds for i in 1:n
        Ct = is_scalar_ct ? c_t : c_t[i]
        yaw_i = is_scalar_yaw ? yaw : yaw[i]
        ti_i  = is_scalar_ti  ? ti  : ti[i]
        ti0_i = is_scalar_ti0 ? ti0 : ti0[i]
        I = sqrt(ti_i^2 + ti0_i^2)
        OPdw = rps[i, 1]
        cosyaw = cos(yaw_i)

        k_y_i = k_a * I + k_b
        k_z_i = k_y_i
        x0_i = x_0[i]

        # Field widths (Eq. 7.2)
        max_term = max(OPdw - x0_i, 0.0)
        min_term = min(OPdw / x0_i, 1.0)
        sig_y_i = max_term * k_y_i + min_term * cosyaw * d_rotor * invsqrt8
        sig_z_i = max_term * k_z_i + min_term * d_rotor * invsqrt8
        sig_y[i] = sig_y_i
        sig_z[i] = sig_z_i

        # Deflection
        Θ = 0.3 * yaw_i / cosyaw * (1 - sqrt(1 - Ct * cosyaw))
        delta_nfw = Θ * min(OPdw, x0_i)
        sqrt_ct = sqrt(Ct)
        delta_fw_1 = Θ / 14.7 * sqrt(cosyaw / (k_y_i * k_z_i * Ct)) * (2.9 + 1.3 * sqrt(1 - Ct) - Ct)
        term = 1.6 * sqrt((8 * sig_y_i * sig_z_i) / (d_rotor^2 * cosyaw))
        num = (1.6 + sqrt_ct) * (term - sqrt_ct)
        den = (1.6 - sqrt_ct) * (term + sqrt_ct)
        arg = num / den
        delta_fw_2 = log(max(arg, eps()))
        blend = 0.5 * sign(OPdw - x0_i) + 0.5
        delta[i, 1] = delta_nfw + blend * delta_fw_1 * delta_fw_2 * d_rotor
        delta[i, 2] = 0.0

        # Potential core
        u_r_0 = (Ct * cosyaw) / (2 * (1 - sqrt(1 - Ct * cosyaw)) * sqrt(1 - Ct))
        ratio_term = max(1 - OPdw / x0_i, 0.0)
        sqrt_u = sqrt(u_r_0)
        pc_y[i] = d_rotor * cosyaw * sqrt_u * ratio_term
        pc_z[i] = d_rotor * sqrt_u * ratio_term
        if OPdw == 0.0
            pc_y[i] = d_rotor * cosyaw
            pc_z[i] = d_rotor
        end
    end

    return nothing
end

"""
    centerline!(deflection::AbstractMatrix,
                states_op, states_t, states_wf, floris::Floris, d_rotor)

Compute the cross-wind wake deflection at the observation points in-place using the Gaussian wake model.

This function calculates the lateral (y) and vertical (z) deflection of the wake centerline 
due to yaw misalignment and other effects. The results are written directly into the provided 
deflection matrix without allocating temporary arrays.

# Output Arguments
- `deflection::AbstractMatrix` (size n×2): Wake deflection components filled in-place
  - Column 1: Lateral deflection Δy `[m]`  
  - Column 2: Vertical deflection Δz `[m]` (always zero in current implementation)

# Input Arguments
- `states_op::AbstractMatrix` (size n×k): Operational point states where n is number of points
  - Column 4 contains downstream distance in wake-aligned coordinates `[m]`
- `states_t::AbstractMatrix`: Turbine state matrix containing:
  - Column 1: Axial induction factor a `[-]`
  - Column 2: Yaw angle `[degrees]`
  - Column 3: Local turbulence intensity TI `[-]`
- `states_wf::AbstractMatrix`: Wind field state matrix containing:
  - Column 3: Ambient turbulence intensity TI₀ `[-]`
- `floris::Floris`: FLORIS model parameters for wake calculations (see [`Floris`](@ref))
- `d_rotor::Real`: Rotor diameter D `[m]`

# Notes
- Only `states_op[:, 4]` (downstream distance) is used from the operational points
- The function internally converts yaw angles from degrees to radians with sign correction
- Thrust coefficient is calculated from axial induction factor using `calcCt`

See also: [`getVars!`](@ref)
"""
function centerline!(deflection::AbstractMatrix,
                     states_op, states_t, states_wf, floris::Floris, d_rotor)
    n = size(states_op, 1)
    # Prepare minimal RPs matrix: only downstream distance (OPdw) is needed (col 1)
    RPs = Matrix{Float64}(undef, n, 3)
    @inbounds begin
        for i in 1:n
            RPs[i, 1] = states_op[i, 4]  # downstream distance in wake coords
            RPs[i, 2] = 0.0
            RPs[i, 3] = 0.0
        end
    end

    # States (vectorised)
    Ct  = calcCt(states_t[:, 1], states_t[:, 2])
    yaw = -deg2rad.(states_t[:, 2])
    TI  = states_t[:, 3]
    TI0 = states_wf[:, 3]

    # Temporary outputs we don't need to keep
    sig_y = Vector{Float64}(undef, n)
    sig_z = Vector{Float64}(undef, n)
    x_0   = Vector{Float64}(undef, n)
    pc_y  = Vector{Float64}(undef, n)
    pc_z  = Vector{Float64}(undef, n)

    # Compute deflection into provided matrix
    getVars!(sig_y, sig_z, x_0, deflection, pc_y, pc_z, RPs, Ct, yaw, TI, TI0, floris, d_rotor)
    return deflection
end

"""
    States

Lightweight container for state variable names and counts used in wind farm simulations.

This mutable struct organizes the state variables into three categories: turbine states,
observation point (OP) states, and wind field (WF) states. It provides both the variable
names and their counts for efficient memory allocation and state management.

# Fields
- `T_names::Vector{String}`: Names of turbine state variables ["a", "yaw", "TI"]
- `Turbine::Int`: Number of turbine state variables (3)
- `OP_names::Vector{String}`: Names of observation point variables ["x0", "y0", "z0", "x1", "y1", "z1"]
- `OP::Int`: Number of observation point variables (6)
- `WF_names::Vector{String}`: Names of wind field variables ["wind_vel", "wind_dir", "TI0"]
- `WF::Int`: Number of wind field variables (3)

# Constructor
Use the default constructor [`States()`](@ref) to create an instance with predefined
variable names and counts currently used by FLORIDyn.jl.

# Example
```julia
states = States()
println("Turbine states: ", states.T_names)
println("Number of turbine states: ", states.Turbine)
```

# State Variable Definitions
## Turbine States (T_names)
- `"a"`: Axial induction factor [-]
- `"yaw"`: Yaw angle [degrees]
- `"TI"`: Local turbulence intensity [-]

## Observation Point States (OP_names)
- `"x0", "y0", "z0"`: Initial position coordinates [m]
- `"x1", "y1", "z1"`: Final position coordinates [m]

## Wind Field States (WF_names)
- `"wind_vel"`: Wind velocity [m/s]
- `"wind_dir"`: Wind direction [degrees]
- `"TI0"`: Ambient turbulence intensity [-]

See also: [`States()`](@ref)
"""
mutable struct States
    T_names::Vector{String}
    Turbine::Int
    OP_names::Vector{String}
    OP::Int
    WF_names::Vector{String}
    WF::Int
end

"""
    States() -> States

Default constructor for the States struct with predefined wind farm simulation variables.

Creates a States instance with commonly used state variable names and their counts
for wind farm simulations using the FLORIS wake model.

# Returns
- `States`: Initialized struct with:
  - Turbine states: ["a", "yaw", "TI"] (count: 3)
  - Observation point states: ["x0", "y0", "z0", "x1", "y1", "z1"] (count: 6)  
  - Wind field states: ["wind_vel", "wind_dir", "TI0"] (count: 3)

# Example
```julia
states = States()
@assert states.Turbine == 3
@assert states.T_names == ["a", "yaw", "TI"]
```

See also: [`States`](@ref)
"""
function States()
    T_names = ["a", "yaw", "TI"]
    Turbine = length(T_names)
    OP_names = ["x0", "y0", "z0", "x1", "y1", "z1"]
    OP = length(OP_names)
    WF_names = ["wind_vel", "wind_dir", "TI0"]
    WF = length(WF_names)
    return States(T_names, Turbine, OP_names, OP, WF_names, WF)
end

"""
    init_states(set::Settings, wf::WindFarm, wind::Wind, init_turb, 
                floris::Floris, sim::Sim) -> Tuple{Matrix, Matrix, Matrix}

Initialize the state arrays for wind farm simulation using the Gaussian wake model.

This function sets up the initial conditions for turbines, observation points, and wind field 
states based on the provided configuration parameters. It computes initial positions, wind 
conditions, and wake properties for each turbine in the wind farm.

# Arguments
- `set::Settings`: Simulation settings containing configuration options for velocity, direction, and turbulence models
- `wf::WindFarm`: Wind farm object containing turbine positions, dimensions, and state arrays (see [`WindFarm`](@ref))
- `wind::Wind`: Wind conditions including velocity, direction, and turbulence intensity data (see [`Wind`](@ref))
- `init_turb`: Initial turbine state parameters (axial induction factor, yaw angle, turbulence intensity)
- `floris::Floris`: FLORIS model parameters for wake calculations (see [`Floris`](@ref))
- `sim::Sim`: Simulation parameters including time step and start time (see [`Sim`](@ref))

# Returns
A tuple `(states_op, states_t, states_wf)` containing:
- `states_op::Matrix`: Observation point states with 3D coordinates and wake positions
- `states_t::Matrix`: Turbine states including control parameters and operational conditions  
- `states_wf::Matrix`: Wind field states with velocity, direction, and turbulence data

# Description
The function performs the following initialization steps for each turbine:
1. Retrieves wind field data (velocity, direction, turbulence intensity) based on the specified input methods
2. Initializes wind field states at all observation points for the turbine
3. Calculates downwind distances for wake coordinate system
4. Sets initial turbine states from provided parameters
5. Computes crosswind wake deflections using the centerline function
6. Transforms coordinates from wake-relative to world coordinate system
7. Updates observation point positions including turbine base and nacelle offsets

# Notes
- Supports multiple wind input methods including interpolation, constant values, and random walk models
- Handles both 3D and 4D wind field configurations (with optional orientation data)
- Uses SOWFA coordinate system conventions for angle transformations
"""
function init_states(set::Settings, wf::WindFarm, wind::Wind, init_turb, floris::Floris, sim::Sim)
    # Unpack state arrays and parameters
    states_op   = copy(wf.States_OP)
    states_t    = copy(wf.States_T)
    states_wf   = copy(wf.States_WF)
    nT          = wf.nT
    nOP         = wf.nOP
    deltaT      = sim.time_step
    startTime   = sim.start_time

    for iT = 1:nT
        # Retrieve wind field data
        if wind.input_vel == "I_and_I"
            u = getWindSpeedT(set.vel_mode, wind.vel, iT, startTime)
        elseif wind.input_vel == "ZOH_wErrorCov"
            u = wind.vel.Init
        elseif wind.input_vel == "RW_with_Mean"
            # For RW_with_Mean during initialization, use a default value
            # The actual RW_with_Mean logic will be handled in the correction phase
            u = 8.0  # Default wind speed for initialization
        else
            u = getWindSpeedT(set.vel_mode, wind.vel, iT, startTime)
        end

        if wind.input_dir == "RW_with_Mean"
            phi_s = wind.dir.Init[iT]
        else
            phi_s = getWindDirT(set.dir_mode, wind.dir, iT, startTime)
        end

        ti = getWindTiT(set.turb_mode, wind.ti, iT, startTime)

        rangeOPs = ((iT-1)*nOP+1):(iT*nOP)

        # Initialize the States of the OPs and turbines
        states_wf[rangeOPs, 1] .= u
        states_wf[rangeOPs, 2] .= phi_s
        states_wf[rangeOPs, 3] .= ti

        # Add orientation if used
        if length(wf.Names_WF) == 4
            states_wf[rangeOPs, 4] .= phi_s
        end

        # Downwind distance (wake coord)
        states_op[rangeOPs, 4] .= (collect(0:(nOP-1)) .* deltaT .* u)

        # Init turbine states
        states_t[rangeOPs, :] .= init_turb[iT, :]'

    # Crosswind position (in-place)
    centerline!(@view(states_op[rangeOPs, 5:6]), states_op[rangeOPs, :], states_t[rangeOPs, :],
            states_wf[rangeOPs, :], floris, wf.D[iT])

        # Convert wind dir in fitting radians
        phi_w = angSOWFA2world.(states_wf[rangeOPs, 2])

        # World coordinate position x0 and y0 including tower base and nacelle pos
        states_op[rangeOPs, 1] .= cos.(phi_w) .* states_op[rangeOPs, 4] .-
                                  sin.(phi_w) .* states_op[rangeOPs, 5] .+   wf.posBase[iT, 1] .+ wf.posNac[iT, 1]
        states_op[rangeOPs, 2] .= sin.(phi_w) .* states_op[rangeOPs, 4] .+
                                  cos.(phi_w) .* states_op[rangeOPs, 5] .+ wf.posBase[iT, 2] .+ wf.posNac[iT, 2]
        states_op[rangeOPs, 3] .= states_op[rangeOPs, 6] .+ wf.posBase[iT, 3] .+ wf.posNac[iT, 3]
    end

    return states_op, states_t, states_wf
end

"""
    getPower(wf::WindFarm, m::AbstractMatrix, floris::Floris, con::Con)

Calculate the power output of wind turbines in a wind farm simulation.

This function computes the power generated by wind turbines based on their operational states,
wind conditions, and control settings. It accounts for yaw angle effects and optional
yaw range constraints using hyperbolic tangent functions for smooth operational limits.

# Arguments
- `wf::WindFarm`: Wind farm object containing turbine states, dimensions, and operational data with current axial induction factors and yaw angles (see [`WindFarm`](@ref))
- `m::Matrix`: Measurement or simulation data matrix where column 3 contains effective wind speeds at turbine locations [m/s]
- `floris::Floris`: FLORIS model parameters containing air density, drivetrain efficiency, and power curve parameters (see [`Floris`](@ref))
- `con::Con`: Controller configuration object with yaw control settings and operational constraints (see [`Con`](@ref))

# Returns
- `P::Vector{Float64}`: Power output for each turbine in the wind farm [W]

# Description
The function calculates power using the standard wind turbine power equation with yaw corrections:

```julia
P = 0.5 × ρ × A × Cp × U³ × η × cos(γ)^p_p × f_yaw_constraints
```

Where:
- `ρ` is air density [`floris.airDen`] [kg/m³]
- `A` is rotor swept area `π × (D/2)²` [m²]
- `Cp` is power coefficient calculated as `4a(1-a)²` [-]
- `U` is effective wind speed from column 3 of matrix `m` [m/s]
- `η` is drivetrain efficiency [`floris.eta`] [-]
- `γ` is yaw angle [rad]
- `p_p` is yaw power exponent [`floris.p_p`] [-]
- `f_yaw_constraints` is optional yaw range constraint factor [-]

The yaw constraint factor is applied when `con.tanh_yaw` is true:
```julia
f_yaw_constraints = [0.5 × tanh((γ_max - γ) × 50) + 0.5] × 
                           [-0.5 × tanh((γ_min - γ) × 50) + 0.5]
```

# Notes
- Power coefficient is calculated from axial induction factor: `Cp = 4a(1-a)²`
- Yaw effects reduce power output according to `cos(γ)^p_p` where `p_p` is typically 1.88
- Optional yaw range constraints use hyperbolic tangent functions with slope factor 50 for smooth transitions
- When `con.tanh_yaw` is enabled, power is smoothly constrained within [`con.yawRangeMin`, `con.yawRangeMax`]
- The constraint functions approach step functions but provide smooth gradients for optimization
- Axial induction factors are extracted from `wf.States_T[wf.StartI, 1]` for current time step
- Yaw angles are converted from degrees to radians internally
"""
function getPower(wf::WindFarm, m::AbstractMatrix, floris::Floris, con::Con)
    a   = wf.States_T[wf.StartI, 1]
    yaw = deg2rad.(wf.States_T[wf.StartI, 2])
    
    Cp = 4a .* (1 .- a).^2
    ueff = m[:, 3]

    if con.tanh_yaw
        P = 0.5 * floris.airDen * (wf.D / 2).^2 * π .* Cp' .* ueff.^3 .* floris.eta .* 
            (cos.(yaw).^floris.p_p)' .* 
            (0.5 * tanh((-yaw + deg2rad.(con.yawRangeMax)) * 50) + 0.5) * 
            (-0.5 * tanh((-yaw + deg2rad.(con.yawRangeMin)) * 50) + 0.5)
    else
        P = 0.5 * floris.airDen * (wf.D / 2).^2 * π .* Cp' .* ueff.^3 .* floris.eta .* 
            (cos.(yaw).^floris.p_p)'
    end

    return P
end

"""
    getUadv(states_op, states_t, states_wf, floris::Floris, d_rotor)

Calculate the advection speed factor based on Zong & Porté-Agel 2020.

The advection speed is defined as:
```
Uadv(x) = 0.5 * (U_inf + U_centre(x))
```

Solved for the ratio Uadv(x)/U_inf:
```
Uadv(x)/Uinf = sqrt(1 - Ct*cos(gamma)/(8*(sig_y*sig_z)/D^2))
```

# Arguments
- `states_op`: Observation point states matrix containing downstream distances in column 4
- `states_t`: Turbine states matrix with axial induction factors (col 1), yaw angles (col 2), and turbulence intensity (col 3)
- `states_wf`: Wind field states matrix with ambient turbulence intensity in column 3
- `floris::Floris`: FLORIS model parameters containing wake expansion coefficients (see [`Floris`](@ref))
- `d_rotor`: Rotor diameter [m]

# Returns
- `Uadv_div_U_inf::Vector{Float64}`: Advection speed factor (ratio of advection speed to freestream wind speed) [-]

# Description
The function calculates the advection speed factor using the following steps:
1. Computes thrust coefficient from axial induction factor
2. Calculates potential core length (x_0) based on turbulence intensity and wake expansion
3. Determines wake field widths (sig_y, sig_z) for Gaussian wake model
4. Applies different formulations inside vs. outside the potential core region
5. Returns the advection speed ratio for wake transport calculations

The advection speed represents the speed at which wake features propagate downstream,
which is typically slower than the freestream wind speed due to velocity deficits.

# References
- Zong, H. and Porté-Agel, F. (2020). A momentum-conserving wake superposition method for wind farm power prediction
"""
function getUadv(states_op, states_t, states_wf, floris::Floris, d_rotor)
    # Parameters
    k_a   = floris.k_a
    k_b   = floris.k_b
    alpha = floris.alpha
    beta  = floris.beta

    # States
    C_T = calcCt(states_t[:, 1], states_t[:, 2])
    yaw = -deg2rad.(states_t[:, 2])
    I = sqrt.(states_t[:, 3].^2 .+ states_wf[:, 3].^2)  # I_f & I_0
    OPdw = states_op[:, 4]

    # Calc x_0 (Core length)
    # [1] Eq. 7.3
    x_0 = (cos.(yaw) .* (1 .+ sqrt.(1 .- C_T)) ./ 
          (sqrt(2) .* (alpha .* I .+ beta .* (1 .- sqrt.(1 .- C_T))))) .* d_rotor

    # Calc k_z and k_y based on I
    # [2] Eq.8
    k_y = k_a .* I .+ k_b
    k_z = k_y

    # Get field width y
    # [1] Eq. 7.2
    # To fit the field width, the value linearly increases from 0 to max for dw
    # positions before x_0
    zs = zeros(size(OPdw))
    sig_y_div_D = max.(OPdw .- x_0, zs) .* k_y ./ d_rotor .+ 
                  min.(OPdw ./ x_0, zs .+ 1) .* cos.(yaw) ./ sqrt(8)

    # Get field width z
    # [1] Eq. 7.2
    sig_z_div_D = max.(OPdw .- x_0, zs) .* k_z ./ d_rotor .+ 
                  min.(OPdw ./ x_0, zs .+ 1) ./ sqrt(8)

    # Calc advection speed ratio
    # U_adv(x) = 0.5*(U_inf + U_cen(x))
    # =>  U_adv(x)/Uinf = 0.5*(1 + U_cen/U_inf)
    
    # Centerspeed outside of the potential core
    # Clamp the argument to prevent negative values under sqrt
    sqrt_arg = max.(1 .- (C_T .* cos.(yaw)) ./ (8 .* sig_y_div_D .* sig_z_div_D), 0.0)
    U_cen_div_U_inf = sqrt.(sqrt_arg)
    
    # Centerspeed inside of the potential core
    inside_core = x_0 .> OPdw
    U_cen_div_U_inf[inside_core] = sqrt.(max.(1 .- C_T[inside_core], 0.0))
    
    Uadv_div_U_inf = 0.5 .* (1 .+ U_cen_div_U_inf)

    return Uadv_div_U_inf
end
