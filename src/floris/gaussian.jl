# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    calcCt(a, _) -> Number or Vector

Calculate the thrust coefficient `ct` for a wind turbine based on the axial induction factor `a`.

# Arguments
- `a`: Axial induction factor, typically between 0 and 0.5 (can be a scalar or vector)
- _: unused parameter

# Returns
- ` ct::Number`: The calculated thrust coefficient.
"""
@inline function calcCt(a, _)
    ct = 4 .* a .* (1 .- a)
    return ct
end

"""
    States

A mutable struct representing the state variables used in the Gaussian wake model. 
This struct is intended to store and update the dynamic properties of the wake during simulation.

# Fields
- `T_names::Vector{String}`: Names of turbine state variables (e.g., "a", "yaw", "TI" for axial induction factor, yaw angle, and turbulence intensity)
- `Turbine::Int`: Number of turbine state variables (length of T_names)
- `OP_names::Vector{String}`: Names of observation point state variables (e.g., "x0", "y0", "z0", "x1", "y1", "z1" for 3D coordinates)
- `OP::Int`: Number of observation point state variables (length of OP_names)
- `WF_names::Vector{String}`: Names of wind field state variables (e.g., "wind_vel", "wind_dir", "TI0" for velocity, direction, and ambient turbulence)
- `WF::Int`: Number of wind field state variables (length of WF_names)

# Description
This struct organizes the naming and counting of different types of state variables used in wind farm simulations:
- Turbine states represent individual turbine properties and control settings
- Observation point states track spatial coordinates for wake calculations
- Wind field states capture environmental conditions affecting the entire wind farm
"""
mutable struct States
    T_names::Vector{String}
    Turbine::Int
    OP_names::Vector{String}
    OP::Int
    WF_names::Vector{String}
    WF::Int
end

function States()
    # Turbine states
    T_names = ["a", "yaw", "TI"]
    Turbine = length(T_names)

    # Observation point states
    OP_names = ["x0", "y0", "z0", "x1", "y1", "z1"]
    OP = length(OP_names)

    # Wind field states
    WF_names = ["wind_vel", "wind_dir", "TI0"]
    WF = length(WF_names)

    return States(T_names, Turbine, OP_names, OP, WF_names, WF)
end

"""
    centerline(states_op, states_t, states_wf, floris, d_rotor) -> Matrix{Float64}

Compute the centerline wake properties for a wind farm simulation.

# Arguments
- `states_op`: Operational states of the turbines (e.g., yaw, pitch, etc.).
- `states_t`: Turbine-specific states (e.g., rotor speed, torque, etc.).
- `states_wf`: Wind farm-level states (e.g., wind direction, wind speed, etc.).
- `floris`: Parameters for the FLORIS wake model.
- `d_rotor`: Rotor diameter or characteristic length scale.

# Returns
- The computed centerline wake properties `delta`, which includes the deflection in the y and z directions.

# Notes
This function is part of the Gaussian wake model implementation.
"""
function centerline(states_op, states_t, states_wf, floris, d_rotor)
    N          = size(states_op, 1)                     # number of rows
    Δ          = Matrix{Float64}(undef, N, 2)        # output
    is8      = 1 / sqrt(8)
    s2         = sqrt(2)
    α, β       = floris.alpha, floris.beta
    k_a, k_b   = floris.k_a, floris.k_b

    for i ∈ 1:N          # ← one tight SIMD loop
        ############# 1. unpack state rows ####################################
        a  = states_t[i, 1]
        b  = states_t[i, 2]
        c  = states_t[i, 3]
        d  = states_wf[i, 3]
        OP = states_op[i, 4]

        ############# 2. basic derived quantities #############################
        C_T     = calcCt(a, b)
        yaw     = -deg2rad(b)
        cosyaw  = cos(yaw)
        I       = sqrt(c*c + d*d)

        ############# 3. core length x₀ #######################################
        x₀ = (cosyaw * (1 + sqrt(1 - C_T)) /
             (s2 * (α * I + β * (1 - sqrt(1 - C_T))))) * d_rotor

        ############# 4. wake expansion (σ_y / σ_z) ###########################
        k_y     = k_a * I + k_b
        k_z     = k_y                    # identical expression

        diff    = OP - x₀
        f1      = max(diff, 0.0)
        f2      = min(OP / x₀, 1.0)

        σ_y = f1 * k_y + f2 * cosyaw * d_rotor * is8
        σ_z = f1 * k_z + f2 * d_rotor    * is8

        ############# 5. deflection angles ####################################
        Θ = 0.3 * yaw / cosyaw * (1 - sqrt(1 - C_T * cosyaw))

        ############# 6. near‑field / far‑field pieces #########################
        Δ_nfw = Θ * min(OP, x₀)

        Δ_fw₁ = Θ / 14.7 * sqrt(cosyaw / (k_y * k_z * C_T)) *
                (2.9 + 1.3 * sqrt(1 - C_T) - C_T)

        num  = (1.6 + sqrt(C_T)) *
               (1.6 * sqrt((8 * σ_y * σ_z) / (d_rotor^2 * cosyaw)) - sqrt(C_T))
        den  = (1.6 - sqrt(C_T)) *
               (1.6 * sqrt((8 * σ_y * σ_z) / (d_rotor^2 * cosyaw)) + sqrt(C_T))

        # The real part of the complex logarithm is log(abs(ratio)), which is what we need.
        Δ_fw₂ = real(log(complex(num / den, 0.0)))

        factor = sign(OP - x₀) / 2 + 0.5        # 0, 0.5, or 1

        Δy = Δ_nfw + factor * Δ_fw₁ * Δ_fw₂ * d_rotor

        ############# 7. store result ##########################################
        Δ[i, 1] = Δy
        Δ[i, 2] = 0.0        # Δz was always zero in the original
    end

    return Δ
end

"""
    centerline!(deflection, states_op, states_t, states_wf, floris, d_rotor)

In-place version of the `centerline` function that calculates the wake centerline deflection.

This function computes the wake centerline deflection for each operational point without 
allocating new memory for the result. The results are stored in the pre-allocated `deflection`
matrix.

# Arguments
- `deflection::AbstractMatrix`: A pre-allocated matrix (`N × 2`) to store the deflection results (Δy, Δz).
- `states_op`, `states_t`, `states_wf`: Matrices containing the operational point, turbine, and wind farm states.
- `floris::Floris`: FLORIS model parameters.
- `d_rotor::Real`: The rotor diameter.

# Returns
- `nothing`: The `deflection` matrix is modified in-place.
"""
function centerline!(deflection::AbstractMatrix, states_op, states_t, states_wf, floris, d_rotor)
    N = size(states_op, 1)
    is8 = 1 / sqrt(8)
    s2 = sqrt(2)
    α, β = floris.alpha, floris.beta
    k_a, k_b = floris.k_a, floris.k_b

    @inbounds for i in 1:N
        # Unpack states
        a = states_t[i, 1]
        b = states_t[i, 2]
        c = states_t[i, 3]
        d = states_wf[i, 3]
        OP = states_op[i, 4]

        # Basic derived quantities
        C_T = calcCt(a, b)
        yaw = -deg2rad(b)
        cosyaw = cos(yaw)
        I = sqrt(c * c + d * d)

        # Core length x₀
        x₀ = (cosyaw * (1 + sqrt(1 - C_T)) /
             (s2 * (α * I + β * (1 - sqrt(1 - C_T))))) * d_rotor

        # Wake expansion (σ_y / σ_z)
        k_y = k_a * I + k_b
        k_z = k_y

        diff = OP - x₀
        f1 = max(diff, 0.0)
        f2 = min(OP / x₀, 1.0)

        σ_y = f1 * k_y + f2 * cosyaw * d_rotor * is8
        σ_z = f1 * k_z + f2 * d_rotor * is8

        # Deflection angles
        Θ = 0.3 * yaw / cosyaw * (1 - sqrt(1 - C_T * cosyaw))

        # Near-field / far-field pieces
        Δ_nfw = Θ * min(OP, x₀)

        Δ_fw₁ = Θ / 14.7 * sqrt(cosyaw / (k_y * k_z * C_T)) *
                (2.9 + 1.3 * sqrt(1 - C_T) - C_T)

        num = (1.6 + sqrt(C_T)) *
              (1.6 * sqrt((8 * σ_y * σ_z) / (d_rotor^2 * cosyaw)) - sqrt(C_T))
        den = (1.6 - sqrt(C_T)) *
              (1.6 * sqrt((8 * σ_y * σ_z) / (d_rotor^2 * cosyaw)) + sqrt(C_T))

        # The real part of the complex logarithm is log(abs(ratio)), which is what we need.
        Δ_fw₂ = real(log(complex(num / den, 0.0)))

        factor = sign(OP - x₀) / 2 + 0.5

        Δy = Δ_nfw + factor * Δ_fw₁ * Δ_fw₂ * d_rotor

        # Store result in the pre-allocated matrix
        deflection[i, 1] = Δy
        deflection[i, 2] = 0.0  # Δz is always zero
    end
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

        # Crosswind position
        states_op[rangeOPs, 5:6] = centerline(states_op[rangeOPs, :], states_t[rangeOPs, :],
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
    getVars(rps::Union{Matrix, Adjoint}, c_t, yaw, ti, ti0, 
            floris::Floris, d_rotor) -> Tuple{Vector, Vector, Vector, Vector, Vector, Vector, Vector}

Compute and return variables related to the Gaussian wake model for wind turbines.

In particular, it calculates the field width, the potential core data and
the deflection. These values are needed for the calculation of the wake shape and speed 
reduction. The values are based of the state of every individual OP.

# Arguments
- `rps`: Matrix of reference points where the variables are evaluated.
- `c_t`: Thrust coefficient(s) for the turbine(s).
- `yaw`: Yaw angle(s) of the turbine(s) in radians or degrees.
- `ti`: Turbulence intensity at the reference points.
- `ti0`: Ambient turbulence intensity.
- `floris::Floris`: FLORIS model parameters containing Gaussian wake model parameters (see [`Floris`](@ref)).
- `d_rotor`: Rotor diameter(s) of the turbine(s).

# Returns
Returns the tuple
- sig_y::Vector: Gaussian variance in y direction (sqrt of)
- sig_z::Vector: Gaussian variance in z direction (sqrt of)
- c_t: Thrust coefficient, same as OP.Ct
- x_0: Potential core length
- delta: Deflection
- pc_y: Potential core boundary in y dir
- pc_z: Potential core boundary in z dir

# SOURCES
- [1] Experimental and theoretical study of wind turbine wakes in yawed conditions - M. Bastankhah and F. Porté-Agel
- [2] Design and analysis of a spatially heterogeneous wake - A. Farrell, J. King et al.
"""
function getVars(rps::Union{Matrix, Adjoint}, c_t, yaw, ti, ti0, floris::Floris, d_rotor)
    # Unpack parameters
    k_a   = floris.k_a
    k_b   = floris.k_b
    alpha = floris.alpha
    beta  = floris.beta

    # States
    I = sqrt.(ti.^2 .+ ti0.^2)
    OPdw = rps[:, 1]

    # Core length x_0
    x_0 = (cos.(yaw) .* (1 .+ sqrt.(1 .- c_t)) ./ 
          (sqrt(2) .* (alpha .* I .+ beta .* (1 .- sqrt.(1 .- c_t))))) .* d_rotor

    # Compute k_y and k_z
    k_y = k_a .* I .+ k_b
    k_z = k_y

    # Helper zero array
    zs = zeros(size(OPdw))

    # sig_y calculation (field width in y)
    sig_y = max.(OPdw .- x_0, zs) .* k_y .+
            min.(OPdw ./ x_0, zs .+ 1) .* cos.(yaw) .* d_rotor / sqrt(8)

    # sig_z calculation (field width in z)
    sig_z = max.(OPdw .- x_0, zs) .* k_z .+
            min.(OPdw ./ x_0, zs .+ 1) .* d_rotor / sqrt(8)

    # Theta
    Theta = 0.3 .* yaw ./ cos.(yaw) .* (1 .- sqrt.(1 .- c_t .* cos.(yaw)))

    # Deflection delta - near wake
    delta_nfw = Theta .* map((opdw, x0) -> min(opdw, x0), OPdw, x_0)

    # delta_fw parts
    delta_fw_1 = Theta ./ 14.7 .* sqrt.(cos.(yaw) ./ (k_y .* k_z .* c_t)) .* 
                 (2.9 .+ 1.3 .* sqrt.(1 .- c_t) .- c_t)

    # Intermediate term
    term = 1.6 .* sqrt.((8 .* sig_y .* sig_z) ./ (d_rotor.^2 .* cos.(yaw)))
    arg = (1.6 .+ sqrt.(c_t)) .* (term .- sqrt.(c_t)) ./ ((1.6 .- sqrt.(c_t)) .* (term .+ sqrt.(c_t)))
    delta_fw_2 = log.(max.(eps(), arg))

    # Condition mask: OPdw > x_0 => 1.0, else 0.0
    mask = (OPdw .> x_0)
    blend = 0.5 .* sign.(OPdw .- x_0) .+ 0.5

    # Total delta in y
    deltaY = delta_nfw .+ blend .* delta_fw_1 .* delta_fw_2 .* d_rotor
    delta = hcat(deltaY, zeros(size(deltaY)))  # [delta_y, delta_z]

    # Potential core
    u_r_0 = (c_t .* cos.(yaw)) ./ 
            (2 .* (1 .- sqrt.(1 .- c_t .* cos.(yaw))) .* sqrt.(1 .- c_t))

    pc_y = d_rotor .* cos.(yaw) .* sqrt.(u_r_0) .* max.(1 .- OPdw ./ x_0, zs)
    pc_z = d_rotor .* sqrt.(u_r_0) .* max.(1 .- OPdw ./ x_0, zs)

    # For points exactly at the rotor plane
    rp = OPdw .== 0
    if sum(rp) > 0
        pc_y[rp] .= d_rotor .* cos.(yaw)
        pc_z[rp] .= d_rotor
    end

    return sig_y, sig_z, c_t, x_0, delta, pc_y, pc_z
end

"""
    runFLORIS(set::Settings, location_t, states_wf, states_t, d_rotor, 
              floris::Floris, windshear::Union{Matrix, WindShear})

Execute the FLORIS (FLOw Redirection and Induction in Steady State) wake model simulation for wind farm analysis.

This function performs a comprehensive wake analysis using the Gaussian wake model to calculate
velocity reductions, turbulence intensity additions, and effective wind speeds at turbine locations.
It accounts for wake interactions, rotor discretization, wind shear effects, and turbulence propagation.

# Arguments
- `set::Settings`: Simulation settings containing configuration options for wind shear modeling
- `location_t`: Matrix of turbine positions [x, y, z] coordinates for each turbine [m]
- `states_wf`: Wind field state matrix containing velocity, direction, and turbulence data
- `states_t`: Turbine state matrix with axial induction factors, yaw angles, and turbulence intensities
- `d_rotor`: Vector of rotor diameters for each turbine [m]
- `floris::Floris`: FLORIS model parameters containing wake model coefficients and rotor discretization settings (see [`Floris`](@ref))
- `windshear`: Wind shear profile data for vertical wind speed variation modeling, either a matrix or of type [`WindShear`](@ref)

# Returns
A tuple `(T_red_arr, T_aTI_arr, T_Ueff, T_weight)` containing:
- `T_red_arr::Vector{Float64}`: Velocity reduction factors for each turbine [-]
- `T_aTI_arr::Union{Vector{Float64}, Nothing}`: Added turbulence intensity from upstream wakes [%]
- `T_Ueff::Union{Float64, Nothing}`: Effective wind speed at the last turbine location [m/s]
- `T_weight::Union{Vector{Float64}, Nothing}`: Gaussian weight factors for wake overlap calculations [-]

# Description
The function performs the following computational steps:

## 1. Rotor Discretization
- Discretizes rotor planes into radial points using the isocell algorithm
- Applies yaw rotation transformations to rotor point coordinates
- Handles both active turbines (d_rotor > 0) and placeholder turbines

## 2. Single Turbine Case
For single turbine simulations, only wind shear effects are calculated without wake interactions.

## 3. Multi-Turbine Wake Analysis
For each upstream turbine affecting downstream turbines:

### Coordinate Transformations
- Transforms rotor points to wake coordinate system aligned with wind direction
- Applies rotational matrices for wind direction and turbine yaw angles
- Filters turbines based on minimum downstream distance (10 rotor diameters)

### Wake Variable Calculations
- Computes wake expansion coefficients, potential core dimensions, and deflection using `getVars`
- Calculates crosswind wake positions and radial distances from wake centerline
- Determines core region boundaries and near/far wake transitions

### Velocity Reduction Modeling
- Applies different deficit models for core region vs. Gaussian wake regions
- Uses velocity deficit superposition for multiple wake interactions
- Accounts for yaw-induced wake deflection and asymmetry

### Turbulence Intensity Enhancement
- Calculates added turbulence intensity using empirical correlations
- Applies Gaussian weighting for spatial distribution of turbulence enhancement
- Uses parameters k_fa, k_fb, k_fc, k_fd for turbulence intensity modeling

## 4. Wind Shear Integration
- Applies vertical wind shear corrections to the downstream turbine
- Uses wind shear profile data for realistic boundary layer effects

# Mathematical Models
The function implements several key wake modeling equations:

**Velocity Deficit**: Based on Gaussian wake theory with yaw corrections
**Deflection**: Uses analytical wake deflection models for yawed turbines  
**Turbulence**: Empirical correlations for wake-added turbulence intensity
**Superposition**: Linear superposition of velocity deficits from multiple wakes

# Notes
- Uses SOWFA (Simulator for Offshore Wind Farm Applications) coordinate conventions
- Implements state-of-the-art Gaussian wake model with yaw considerations
- Supports both research and engineering applications for wind farm optimization
- Computational complexity scales as O(N²) for N turbines due to wake interactions
- Requires proper initialization of turbine states and wind field conditions

# References
- Bastankhah, M. and Porté-Agel, F. (2016). Experimental and theoretical study of wind turbine wakes in yawed conditions
- Niayifar, A. and Porté-Agel, F. (2016). Analytical modeling of wind farms: A new approach for power prediction
"""
function runFLORIS(set::Settings, location_t, states_wf, states_t, d_rotor, floris::Floris, 
                   windshear::Union{Matrix, WindShear}; alloc=nothing)
    a = @allocated if d_rotor[end] > 0
        RPl, RPw = discretizeRotor(floris.rotor_points)
    else
        RPl = [0.0 0.0 0.0]
        RPw = [1.0]
    end
    if ! isnothing(alloc)
        alloc.if1 += a
        alloc.n += 1
    end
    # Yaw rotation for last turbine
    tmp_yaw = deg2rad(states_t[end, 2])
    R = SA[cos(tmp_yaw)  sin(tmp_yaw)  0.0;
          -sin(tmp_yaw)  cos(tmp_yaw)  0.0;
           0.0           0.0           1.0]

    a = @allocated RPl = (R * (RPl .* d_rotor[end])')' .+ location_t[end, :]'
    if ! isnothing(alloc)
        alloc.expr1 += a
    end

    a = @allocated if length(d_rotor) == 1
        redShear = getWindShearT(set.shear_mode, windshear, RPl[:, 3] ./ location_t[3])
        T_red_arr = RPw' * redShear
        T_aTI_arr, T_Ueff, T_weight = nothing, nothing, nothing
        return T_red_arr, T_aTI_arr, T_Ueff, T_weight
    end
    if ! isnothing(alloc)
        alloc.if2 += a
    end

    # Initialize outputs
    a = @allocated begin
        nT = length(d_rotor)
        T_red_arr = ones(nT)
        T_aTI_arr = zeros(nT - 1)
        T_weight = zeros(nT - 1)
    end
    if ! isnothing(alloc)
        alloc.expr2 += a
    end

    a = @allocated for iT in 1:(nT - 1)

        tmp_phi = size(states_wf,2) == 4 ? angSOWFA2world(states_wf[iT, 4]) :
                                           angSOWFA2world(states_wf[iT, 2])

        tmp_RPs = RPl .- location_t[iT, :]'
        R_phi = [cos(tmp_phi)  sin(tmp_phi)  0.0;
                -sin(tmp_phi)  cos(tmp_phi)  0.0;
                 0.0           0.0           1.0]

        tmp_RPs = (R_phi * tmp_RPs')'

        if tmp_RPs[1, 1] <= 10
            continue
        end

        a = states_t[iT, 1]
        yaw_deg = states_t[iT, 2]
        yaw = -deg2rad(yaw_deg)
        TI = states_t[iT, 3]
        Ct = calcCt(a, yaw_deg)
        TI0 = states_wf[iT, 3]

        sig_y, sig_z, C_T, x_0, delta, pc_y, pc_z = getVars(
            tmp_RPs, Ct, yaw, TI, TI0, floris, d_rotor[iT]
        )

        cw_y = tmp_RPs[:, 2] .- delta[:, 1]
        cw_z = tmp_RPs[:, 3] .- delta[:, 2]
        phi_cw = atan.(cw_z, cw_y)
        r_cw = sqrt.(cw_y.^2 .+ cw_z.^2)

        core = (r_cw .< sqrt.((0.5 .* pc_y .* cos.(phi_cw)).^2 .+
                              (0.5 .* pc_z .* sin.(phi_cw)).^2)) .|
               (tmp_RPs[:, 1] .== 0.0)

        nw = tmp_RPs[:, 1] .< x_0

        tmp_RPs_r = zeros(size(RPw))
        tmp_RPs_r[core] .= 1 .- sqrt(1 .- C_T)

        fw = .!nw
        gaussAbs = zeros(size(RPw))
        gaussAbs[nw] .= 1 .- sqrt(1 .- C_T)
        gaussAbs[fw] .= 1 .- sqrt.(1 .- C_T .* cos(yaw) ./ (8 .* sig_y[fw] .* sig_z[fw] ./ d_rotor[iT]^2))

        gaussWght = ones(size(RPw))
        not_core = .!core
        if any(not_core)
            exp_y = @. exp(-0.5 * ((cw_y[not_core] - cos(phi_cw[not_core]) .* pc_y[not_core] * 0.5) ./ sig_y[not_core])^2)
            exp_z = @. exp(-0.5 * ((cw_z[not_core] - sin(phi_cw[not_core]) .* pc_z[not_core] * 0.5) ./ sig_z[not_core])^2)

            gaussWght[not_core] .= exp_y .* exp_z
            tmp_RPs_r[not_core] .= gaussAbs[not_core] .* gaussWght[not_core]
        end

        T_weight[iT] = sum(gaussWght)
        T_red_arr[iT] = 1 .- dot(RPw, tmp_RPs_r)

        # Added TI
        T_addedTI_tmp = floris.k_fa * (
            a^floris.k_fb *
            TI0^floris.k_fc *
            (mean(tmp_RPs[:, 1]) / d_rotor[iT])^floris.k_fd
        )

        TIexp = floris.TIexp
        exp_y = @. exp(-0.5 * ((cw_y - cos(phi_cw) .* pc_y * 0.5) ./ (TIexp .* sig_y))^2)
        exp_z = @. exp(-0.5 * ((cw_z - sin(phi_cw) .* pc_z * 0.5) ./ (TIexp .* sig_z))^2)

        T_aTI_arr[iT] = T_addedTI_tmp * dot(RPw, exp_y .* exp_z)
    end
    if ! isnothing(alloc)
        alloc.for1 += a
    end

    a = @allocated begin
        redShear = getWindShearT(set.shear_mode, windshear, RPl[:, 3] ./ location_t[end, 3])
        T_red_arr[end] = dot(RPw, redShear)

        T_red = prod(T_red_arr)
        T_Ueff = states_wf[end, 1] * T_red
    end

    if ! isnothing(alloc)
        alloc.expr3 += a
    end

    return T_red_arr, T_aTI_arr, T_Ueff, T_weight
end

"""
In-place, non-allocating version of runFLORIS for performance-critical applications.

This function performs the same wake analysis as runFLORIS but uses preallocated buffers
and in-place operations to minimize memory allocations during simulation loops.

# Arguments
- `T_red_arr`: Preallocated output array for velocity reduction factors
- `T_aTI_arr`: Preallocated output array for added turbulence intensity
- `T_weight`: Preallocated output array for Gaussian weights
- `buffers`: Struct containing all preallocated temporary arrays
- `set::Settings`: Simulation settings
- `location_t`: Matrix of turbine positions [x, y, z] coordinates
- `states_wf`: Wind field state matrix
- `states_t`: Turbine state matrix
- `d_rotor`: Vector of rotor diameters
- `floris::Floris`: FLORIS model parameters
- `windshear`: Wind shear profile data

# Returns
- `T_Ueff::Float64`: Effective wind speed at the last turbine

All outputs are written to the provided preallocated arrays in-place.
"""
function runFLORIS!(T_red_arr::Vector{Float64}, T_aTI_arr::Vector{Float64}, T_weight::Vector{Float64}, 
                    buffers, set::Settings, location_t, states_wf, states_t, d_rotor, 
                    floris::Floris, windshear::Union{Matrix, WindShear})
    
    nT = length(d_rotor)
    
    # Initialize output arrays in-place
    fill!(T_red_arr, 1.0)
    fill!(T_aTI_arr, 0.0)
    fill!(T_weight, 0.0)
    
    # Get rotor discretization (potentially allocating, but only once)
    if d_rotor[end] > 0
        RPl, RPw = discretizeRotor(floris.rotor_points)
    else
        # Use single point for inactive turbine
        copyto!(buffers.RPl_single, [0.0 0.0 0.0])
        buffers.RPw_single[1] = 1.0
        RPl, RPw = view(buffers.RPl_single, 1:1, :), view(buffers.RPw_single, 1:1)
    end
    
    # Yaw rotation for last turbine (in-place using buffers)
    tmp_yaw = deg2rad(states_t[end, 2])
    cos_yaw, sin_yaw = cos(tmp_yaw), sin(tmp_yaw)
    
    # Rotation matrix computation in-place
    buffers.R[1,1] = cos_yaw;  buffers.R[1,2] = sin_yaw;  buffers.R[1,3] = 0.0
    buffers.R[2,1] = -sin_yaw; buffers.R[2,2] = cos_yaw;  buffers.R[2,3] = 0.0
    buffers.R[3,1] = 0.0;      buffers.R[3,2] = 0.0;      buffers.R[3,3] = 1.0
    
    # Apply rotation and scaling in-place
    n_points = size(RPl, 1)
    for i in 1:n_points
        # Scale by rotor diameter
        for j in 1:3
            buffers.tmp_vec[j] = RPl[i, j] * d_rotor[end]
        end
        
        # Apply rotation
        for j in 1:3
            buffers.RPl_transformed[i, j] = 0.0
            for k in 1:3
                buffers.RPl_transformed[i, j] += buffers.R[j, k] * buffers.tmp_vec[k]
            end
            # Add location offset
            buffers.RPl_transformed[i, j] += location_t[end, j]
        end
    end
    
    # Single turbine case
    if nT == 1
        # Compute wind shear in-place
        for i in 1:n_points
            buffers.height_ratios[i] = buffers.RPl_transformed[i, 3] / location_t[1, 3]
        end
        buffers.redShear_view = view(buffers.redShear, 1:n_points)
        getWindShearT!(buffers.redShear_view, set.shear_mode, windshear, view(buffers.height_ratios, 1:n_points))
        
        T_red_arr[1] = dot(view(RPw, 1:n_points), buffers.redShear_view)
        return states_wf[end, 1] * T_red_arr[1]
    end
    
    # Multi-turbine wake analysis
    for iT in 1:(nT - 1)
        # Get wind direction
        tmp_phi = if size(states_wf, 2) == 4
            angSOWFA2world(states_wf[iT, 4])
        else
            angSOWFA2world(states_wf[iT, 2])
        end
        
        # Compute relative positions in-place
        for i in 1:n_points
            for j in 1:3
                buffers.tmp_RPs[i, j] = buffers.RPl_transformed[i, j] - location_t[iT, j]
            end
        end
        
        # Apply wind direction rotation in-place
        cos_phi, sin_phi = cos(tmp_phi), sin(tmp_phi)
        buffers.R_phi[1,1] = cos_phi;  buffers.R_phi[1,2] = sin_phi;  buffers.R_phi[1,3] = 0.0
        buffers.R_phi[2,1] = -sin_phi; buffers.R_phi[2,2] = cos_phi;  buffers.R_phi[2,3] = 0.0
        buffers.R_phi[3,1] = 0.0;      buffers.R_phi[3,2] = 0.0;      buffers.R_phi[3,3] = 1.0
        
        for i in 1:n_points
            for j in 1:3
                buffers.tmp_vec[j] = buffers.tmp_RPs[i, j]
            end
            for j in 1:3
                buffers.tmp_RPs[i, j] = 0.0
                for k in 1:3
                    buffers.tmp_RPs[i, j] += buffers.R_phi[j, k] * buffers.tmp_vec[k]
                end
            end
        end
        
        # Skip if too close (first point check for efficiency)
        if buffers.tmp_RPs[1, 1] <= 10.0
            continue
        end
        
        # Extract turbine states
        a = states_t[iT, 1]
        yaw_deg = states_t[iT, 2]
        yaw = -deg2rad(yaw_deg)
        TI = states_t[iT, 3]
        Ct = calcCt(a, yaw_deg)
        TI0 = states_wf[iT, 3]
        
        # Get wake variables in-place using buffers
        getVars!(buffers.sig_y, buffers.sig_z, buffers.delta, buffers.pc_y, buffers.pc_z,
                 view(buffers.tmp_RPs, 1:n_points, :), Ct, yaw, TI, TI0, floris, d_rotor[iT])
        
        # Compute crosswind positions in-place
        for i in 1:n_points
            buffers.cw_y[i] = buffers.tmp_RPs[i, 2] - buffers.delta[i, 1]
            buffers.cw_z[i] = buffers.tmp_RPs[i, 3] - buffers.delta[i, 2]
            buffers.phi_cw[i] = atan(buffers.cw_z[i], buffers.cw_y[i])
            buffers.r_cw[i] = sqrt(buffers.cw_y[i]^2 + buffers.cw_z[i]^2)
        end
        
        # Compute core region and wake effects in-place
        computeWakeEffects!(view(buffers.tmp_RPs_r, 1:n_points), view(buffers.gaussWght, 1:n_points),
                           view(buffers.cw_y, 1:n_points), view(buffers.cw_z, 1:n_points), 
                           view(buffers.phi_cw, 1:n_points), view(buffers.r_cw, 1:n_points),
                           view(buffers.tmp_RPs, 1:n_points, :), view(buffers.pc_y, 1:n_points), 
                           view(buffers.pc_z, 1:n_points), view(buffers.sig_y, 1:n_points), 
                           view(buffers.sig_z, 1:n_points), Ct, yaw, d_rotor[iT], buffers.x_0[1])
        
        T_weight[iT] = sum(view(buffers.gaussWght, 1:n_points))
        T_red_arr[iT] = 1.0 - dot(view(RPw, 1:n_points), view(buffers.tmp_RPs_r, 1:n_points))
        
        # Added turbulence intensity computation in-place
        T_addedTI_tmp = floris.k_fa * (a^floris.k_fb * TI0^floris.k_fc * 
                                      (sum(view(buffers.tmp_RPs, 1:n_points, 1)) / n_points / d_rotor[iT])^floris.k_fd)
        
        computeAddedTI!(view(buffers.exp_y, 1:n_points), view(buffers.exp_z, 1:n_points),
                       view(buffers.cw_y, 1:n_points), view(buffers.cw_z, 1:n_points),
                       view(buffers.phi_cw, 1:n_points), view(buffers.pc_y, 1:n_points),
                       view(buffers.pc_z, 1:n_points), view(buffers.sig_y, 1:n_points),
                       view(buffers.sig_z, 1:n_points), floris.TIexp)
        
        T_aTI_arr[iT] = T_addedTI_tmp * dot(view(RPw, 1:n_points), 
                                          view(buffers.exp_y, 1:n_points) .* view(buffers.exp_z, 1:n_points))
    end
    
    # Final wind shear computation for last turbine
    for i in 1:n_points
        buffers.height_ratios[i] = buffers.RPl_transformed[i, 3] / location_t[end, 3]
    end
    buffers.redShear_view = view(buffers.redShear, 1:n_points)
    getWindShearT!(buffers.redShear_view, set.shear_mode, windshear, view(buffers.height_ratios, 1:n_points))
    
    T_red_arr[end] = dot(view(RPw, 1:n_points), buffers.redShear_view)
    
    # Compute effective wind speed
    T_red = 1.0
    for i in 1:nT
        T_red *= T_red_arr[i]
    end
    T_Ueff = states_wf[end, 1] * T_red
    
    return T_Ueff
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

"""
Buffer struct for non-allocating runFLORIS! function.

Contains all preallocated arrays needed for the wake calculations.
All arrays should be sized according to the maximum expected rotor discretization points.
"""
mutable struct FLORISBuffers
    # Rotation matrices
    R::Matrix{Float64}              # 3x3 rotation matrix for yaw
    R_phi::Matrix{Float64}          # 3x3 rotation matrix for wind direction
    
    # Rotor point arrays
    RPl_single::Matrix{Float64}     # 1x3 single point for inactive turbines
    RPw_single::Vector{Float64}     # 1-element weight array
    RPl_transformed::Matrix{Float64} # Transformed rotor points
    tmp_RPs::Matrix{Float64}        # Temporary rotor points array
    
    # Temporary vectors
    tmp_vec::Vector{Float64}        # 3-element temporary vector
    height_ratios::Vector{Float64}  # Height ratio calculations
    redShear::Vector{Float64}       # Wind shear reduction factors
    redShear_view::SubArray{Float64,1,Vector{Float64},Tuple{UnitRange{Int64}},true}
    
    # Wake calculation arrays
    sig_y::Vector{Float64}          # Wake width in y direction
    sig_z::Vector{Float64}          # Wake width in z direction
    delta::Matrix{Float64}          # Wake deflection (nx2)
    pc_y::Vector{Float64}           # Potential core boundary y
    pc_z::Vector{Float64}           # Potential core boundary z
    x_0::Vector{Float64}            # Core length (1-element for scalar storage)
    
    # Crosswind position arrays
    cw_y::Vector{Float64}           # Crosswind y positions
    cw_z::Vector{Float64}           # Crosswind z positions
    phi_cw::Vector{Float64}         # Crosswind angles
    r_cw::Vector{Float64}           # Crosswind radial distances
    
    # Wake effect arrays
    tmp_RPs_r::Vector{Float64}      # Wake reduction factors
    gaussWght::Vector{Float64}      # Gaussian weights
    
    # Turbulence intensity arrays
    exp_y::Vector{Float64}          # Exponential terms in y
    exp_z::Vector{Float64}          # Exponential terms in z
end

"""
Create FLORISBuffers with arrays sized for maximum rotor discretization points.
"""
function FLORISBuffers(max_points::Int)
    return FLORISBuffers(
        zeros(3, 3),                    # R
        zeros(3, 3),                    # R_phi
        zeros(1, 3),                    # RPl_single
        ones(1),                        # RPw_single
        zeros(max_points, 3),           # RPl_transformed
        zeros(max_points, 3),           # tmp_RPs
        zeros(3),                       # tmp_vec
        zeros(max_points),              # height_ratios
        zeros(max_points),              # redShear
        view(zeros(max_points), 1:1),   # redShear_view (placeholder)
        zeros(max_points),              # sig_y
        zeros(max_points),              # sig_z
        zeros(max_points, 2),           # delta
        zeros(max_points),              # pc_y
        zeros(max_points),              # pc_z
        zeros(1),                       # x_0
        zeros(max_points),              # cw_y
        zeros(max_points),              # cw_z
        zeros(max_points),              # phi_cw
        zeros(max_points),              # r_cw
        zeros(max_points),              # tmp_RPs_r
        zeros(max_points),              # gaussWght
        zeros(max_points),              # exp_y
        zeros(max_points)               # exp_z
    )
end

"""
In-place version of getVars function that writes results to preallocated buffers.
"""
function getVars!(sig_y, sig_z, delta, pc_y, pc_z, rps, c_t, yaw, ti, ti0, floris, d_rotor)
    # Unpack parameters
    k_a = floris.k_a
    k_b = floris.k_b
    alpha = floris.alpha
    beta = floris.beta
    
    # States
    I = sqrt(ti^2 + ti0^2)
    
    # Core length x_0
    x_0 = (cos(yaw) * (1 + sqrt(1 - c_t)) / 
          (sqrt(2) * (alpha * I + beta * (1 - sqrt(1 - c_t))))) * d_rotor
    
    # Compute k_y and k_z
    k_y = k_a * I + k_b
    k_z = k_y
    
    # Process each point
    n_points = size(rps, 1)
    for i in 1:n_points
        OPdw = rps[i, 1]
        
        # sig_y calculation (field width in y)
        sig_y[i] = max(OPdw - x_0, 0.0) * k_y +
                   min(OPdw / x_0, 1.0) * cos(yaw) * d_rotor / sqrt(8)
        
        # sig_z calculation (field width in z)
        sig_z[i] = max(OPdw - x_0, 0.0) * k_z +
                   min(OPdw / x_0, 1.0) * d_rotor / sqrt(8)
    end
    
    # Theta
    Theta = 0.3 * yaw / cos(yaw) * (1 - sqrt(1 - c_t * cos(yaw)))
    
    # Deflection delta calculation
    for i in 1:n_points
        OPdw = rps[i, 1]
        
        # Near wake deflection
        delta_nfw = Theta * min(OPdw, x_0)
        
        # Far wake deflection parts
        delta_fw_1 = Theta / 14.7 * sqrt(cos(yaw) / (k_y * k_z * c_t)) * 
                     (2.9 + 1.3 * sqrt(1 - c_t) - c_t)
        
        # Intermediate term
        term = 1.6 * sqrt((8 * sig_y[i] * sig_z[i]) / (d_rotor^2 * cos(yaw)))
        arg = (1.6 + sqrt(c_t)) * (term - sqrt(c_t)) / ((1.6 - sqrt(c_t)) * (term + sqrt(c_t)))
        delta_fw_2 = log(max(eps(), arg))
        
        # Blending factor
        blend = 0.5 * sign(OPdw - x_0) + 0.5
        
        # Total delta in y
        delta[i, 1] = delta_nfw + blend * delta_fw_1 * delta_fw_2 * d_rotor
        delta[i, 2] = 0.0  # delta_z is always zero
    end
    
    # Potential core
    u_r_0 = (c_t * cos(yaw)) / (2 * (1 - sqrt(1 - c_t * cos(yaw))) * sqrt(1 - c_t))
    
    for i in 1:n_points
        OPdw = rps[i, 1]
        pc_y[i] = d_rotor * cos(yaw) * sqrt(u_r_0) * max(1 - OPdw / x_0, 0.0)
        pc_z[i] = d_rotor * sqrt(u_r_0) * max(1 - OPdw / x_0, 0.0)
        
        # For points exactly at the rotor plane
        if OPdw == 0.0
            pc_y[i] = d_rotor * cos(yaw)
            pc_z[i] = d_rotor
        end
    end
end

"""
In-place computation of wake effects including core region detection and Gaussian weighting.
"""
function computeWakeEffects!(tmp_RPs_r, gaussWght, cw_y, cw_z, phi_cw, r_cw, tmp_RPs, 
                            pc_y, pc_z, sig_y, sig_z, C_T, yaw, d_rotor, x_0)
    n_points = length(tmp_RPs_r)
    
    # Initialize arrays
    fill!(tmp_RPs_r, 0.0)
    fill!(gaussWght, 1.0)
    
    for i in 1:n_points
        OPdw = tmp_RPs[i, 1]
        
        # Check if in core region
        core_boundary = sqrt((0.5 * pc_y[i] * cos(phi_cw[i]))^2 + 
                            (0.5 * pc_z[i] * sin(phi_cw[i]))^2)
        is_core = (r_cw[i] < core_boundary) || (OPdw == 0.0)
        is_near_wake = OPdw < x_0
        
        if is_core
            tmp_RPs_r[i] = 1 - sqrt(1 - C_T)
        else
            # Gaussian wake region
            if is_near_wake
                gauss_abs = 1 - sqrt(1 - C_T)
            else
                gauss_abs = 1 - sqrt(1 - C_T * cos(yaw) / (8 * sig_y[i] * sig_z[i] / d_rotor^2))
            end
            
            # Gaussian weights
            exp_y = exp(-0.5 * ((cw_y[i] - cos(phi_cw[i]) * pc_y[i] * 0.5) / sig_y[i])^2)
            exp_z = exp(-0.5 * ((cw_z[i] - sin(phi_cw[i]) * pc_z[i] * 0.5) / sig_z[i])^2)
            
            gaussWght[i] = exp_y * exp_z
            tmp_RPs_r[i] = gauss_abs * gaussWght[i]
        end
    end
end

"""
In-place computation of added turbulence intensity exponential terms.
"""
function computeAddedTI!(exp_y, exp_z, cw_y, cw_z, phi_cw, pc_y, pc_z, sig_y, sig_z, TIexp)
    n_points = length(exp_y)
    
    for i in 1:n_points
        exp_y[i] = exp(-0.5 * ((cw_y[i] - cos(phi_cw[i]) * pc_y[i] * 0.5) / (TIexp * sig_y[i]))^2)
        exp_z[i] = exp(-0.5 * ((cw_z[i] - sin(phi_cw[i]) * pc_z[i] * 0.5) / (TIexp * sig_z[i]))^2)
    end
end

"""
In-place version of getWindShearT function that writes results to preallocated output array.

# Arguments
- `output`: Preallocated array to store wind shear reduction factors
- `shear_mode`: Shear model type
- `windshear`: Wind shear parameters or data
- `height_ratios`: Normalized height ratios (h/href)

# Returns
Nothing, results are written to the output array in-place.
"""
function getWindShearT!(output, shear_mode, windshear, height_ratios)
    n = length(height_ratios)
    
    if shear_mode == Shear_PowerLaw()
        alpha = windshear.alpha
        for i in 1:n
            output[i] = height_ratios[i]^alpha
        end
    elseif shear_mode == Shear_Interpolation()
        # For interpolation mode, just call the allocating version for each point
        # This could be optimized further if needed
        for i in 1:n
            output[i] = getWindShearT(shear_mode, windshear, height_ratios[i])
        end
    else
        # Default: no wind shear
        fill!(output, 1.0)
    end
    
    return nothing
end

