# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    calcCt(a, _)

Calculate the thrust coefficient (Ct) for a wind turbine based on the axial induction factor `a`.

# Arguments
- `a::Number`: Axial induction factor, typically between 0 and 0.5.
- _: unused parameter

# Returns
- `Ct::Number`: The calculated thrust coefficient.
"""
function calcCt(a, _)
    Ct = 4 .* a .* (1 .- a)
    return Ct
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
    centerline(states_op, states_t, states_wf, floris, d_rotor)

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
This function is part of the Gaussian wake model implementation for wind farm simulations using the FLORIDyn.jl package.
"""
function centerline(states_op, states_t, states_wf, floris, d_rotor)
    # Parameters
    k_a   = floris.k_a
    k_b   = floris.k_b
    alpha = floris.alpha
    beta  = floris.beta

    # States
    C_T   = calcCt(states_t[:,1], states_t[:,2])
    yaw   = .-deg2rad.(states_t[:,2])
    I     = sqrt.(states_t[:,3].^2 .+ states_wf[:,3].^2)
    OPdw  = states_op[:,4]

    # Calc x_0 (Core length)
    x_0 = (cos.(yaw) .* (1 .+ sqrt.(1 .- C_T)) ./ 
         (sqrt(2) .* (alpha .* I .+ beta .* (1 .- sqrt.(1 .- C_T))))) .* d_rotor

    # Calc k_z and k_y based on I
    k_y = k_a .* I .+ k_b
    k_z = k_y

    # Get field width y
    zs = zeros(size(OPdw))
    sig_y = max.(OPdw .- x_0, zs) .* k_y .+
        min.(OPdw ./ x_0, zs .+ 1) .* cos.(yaw) .* d_rotor ./ sqrt(8)

    # Get field width z
    sig_z = max.(OPdw .- x_0, zs) .* k_z .+
        min.(OPdw ./ x_0, zs .+ 1) .* d_rotor ./ sqrt(8)

    # Calc Theta
    Theta = 0.3 .* yaw ./ cos.(yaw) .* (1 .- sqrt.(1 .- C_T .* cos.(yaw)))

    # Calc Delta/Deflection
    delta_nfw = Theta .* min.(OPdw, x_0)

    delta_fw_1 = Theta ./ 14.7 .* sqrt.(cos.(yaw) ./ (k_y .* k_z .* C_T)) .* (2.9 .+ 1.3 .* sqrt.(1 .- C_T) .- C_T)
    delta_fw_2 = log.(Complex.(
        (1.6 .+ sqrt.(C_T)) .* 
        (1.6 .* sqrt.((8 .* sig_y .* sig_z) ./ (d_rotor^2 .* cos.(yaw))) .- sqrt.(C_T)) ./
        ((1.6 .- sqrt.(C_T)) .*
        (1.6 .* sqrt.((8 .* sig_y .* sig_z) ./ (d_rotor^2 .* cos.(yaw))) .+ sqrt.(C_T)))
    ))
    # println("delta_fw_2: ", delta_fw_2)
    # Use signbit and broadcasting for the sign/(2) + 0.5 logic
    factor = (sign.(OPdw .- x_0) ./ 2 .+ 0.5)

    deltaY = delta_nfw .+ factor .* delta_fw_1 .* delta_fw_2 .* d_rotor
    
    # Deflection in y and z direction
    delta = hcat(deltaY, zeros(size(deltaY)))
    
    return delta
end

"""
    init_states(set::Settings, wf::WindFarm, wind::Wind, init_turb, floris::Floris, sim::Sim)

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
        elseif wind.input_vel in ["ZOH_wErrorCov", "RW_with_Mean"]
            u = wind.vel.Init
        else
            u = getWindSpeedT(set.vel_mode, wind.vel, iT, startTime)
        end

        if wind.input_dir == "RW_with_Mean"
            phi_s = wind.dir.Init
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
        states_t[rangeOPs, :] = ones(nOP, 1) * init_turb[iT, :]'

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
    getVars(rps::Union{Matrix, Adjoint}, c_t, yaw, ti, ti0, floris::Floris, d_rotor)

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

    delta_fw_2 = log.(((1.6 .+ sqrt.(c_t)) .* (term .- sqrt.(c_t))) ./ 
                      ((1.6 .- sqrt.(c_t)) .* (term .+ sqrt.(c_t))))

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
    runFLORIS(set::Settings, location_t, states_wf, states_t, d_rotor, floris::Floris, windshear)

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
- `windshear::Matrix`: Wind shear profile data for vertical wind speed variation modeling

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
function runFLORIS(set::Settings, location_t, states_wf, states_t, d_rotor, floris::Floris, windshear)
    # Main.@infiltrate
    if d_rotor[end] > 0
        RPl, RPw = discretizeRotor(floris.rotor_points)
    else
        RPl = [0.0 0.0 0.0]
        RPw = [1.0]
    end
    # Yaw rotation for last turbine
    tmp_yaw = deg2rad(states_t[end, 2])
    R = [cos(tmp_yaw)  sin(tmp_yaw)  0.0;
        -sin(tmp_yaw)  cos(tmp_yaw)  0.0;
         0.0           0.0           1.0]

    RPl = (R * (RPl .* d_rotor[end])')' .+ location_t[end, :]'

    if length(d_rotor) == 1
        redShear = getWindShearT(set.shear_mode, windshear, RPl[:, 3] ./ location_t[3])
        T_red_arr = RPw' * redShear
        T_aTI_arr, T_Ueff, T_weight = nothing, nothing, nothing
        return T_red_arr, T_aTI_arr, T_Ueff, T_weight
    end

    # Initialize outputs
    nT = length(d_rotor)
    T_red_arr = ones(nT)
    T_aTI_arr = zeros(nT - 1)
    T_weight = zeros(nT - 1)

    for iT in 1:(nT - 1)

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

    redShear = getWindShearT(set.shear_mode, windshear, RPl[:, 3] ./ location_t[end, 3])
    T_red_arr[end] = dot(RPw, redShear)

    T_red = prod(T_red_arr)
    T_Ueff = states_wf[end, 1] * T_red

    return T_red_arr, T_aTI_arr, T_Ueff, T_weight
end

"""
    getPower(wf::WindFarm, m::Matrix, floris::Floris, con::Con)

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

```
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
```
f_yaw_constraints = [0.5 × tanh((γ_max - γ) × 50) + 0.5] × [-0.5 × tanh((γ_min - γ) × 50) + 0.5]
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
function getPower(wf::WindFarm, m::Matrix, floris::Floris, con::Con)
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


