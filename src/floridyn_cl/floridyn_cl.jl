# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    angSOWFA2world(deg_SOWFA)

Convert wind direction angle from SOWFA convention to world coordinate system.

This function performs coordinate transformation between different wind direction conventions
used in wind farm simulations. SOWFA (Simulator fOr Wind Farm Applications) uses a different
angular reference system than the standard world coordinate system used in calculations.

# Arguments
- `deg_SOWFA::Real`: Wind direction angle in SOWFA convention [degrees]

# Returns
- `rad_World::Float64`: Wind direction angle in world coordinate system [radians]

# Coordinate System Conversion
The transformation follows the relationship:
```
θ_world = 270° - θ_SOWFA
```

## SOWFA Convention
- Wind direction angles are defined clockwise from a reference direction

## World Convention  
- Wind direction angles are defined counterclockwise for mathematical calculations
- Standard convention used in wake models and analytical computations

# Mathematical Description
The conversion process:
1. **Angular transformation**: `deg_World = 270 - deg_SOWFA`
2. **Unit conversion**: `rad_World = deg2rad(deg_World)`

The 270° offset accounts for the difference between clockwise (SOWFA) and 
counterclockwise (world) angular conventions.

# Examples
```julia
# Convert 90° SOWFA direction to world coordinates
world_angle = angSOWFA2world(90.0)  # Returns 3.141592... (180° in radians)

# Convert 0° SOWFA direction  
world_angle = angSOWFA2world(0.0)   # Returns 4.712388... (270° in radians)
```

# Notes
- The function handles the sign convention difference between coordinate systems
- Output is always in radians for use in trigonometric calculations
- This transformation is essential for proper wake modeling in wind farm simulations
- The 270° offset ensures proper alignment between SOWFA and mathematical conventions
"""
@inline function angSOWFA2world(deg_SOWFA)
    deg_World = 270.0 - deg_SOWFA
    rad_World = deg2rad(deg_World)
    return rad_World
end

"""
    initSimulation(wf::Union{Nothing, WindFarm}, sim::Sim)

Initialize or load a wind farm simulation state based on simulation settings.

This function handles the initialization phase of a wind farm simulation by either saving 
the current initialized state to disk or loading a previously saved state, depending on 
the simulation configuration.

# Arguments
- `wf::Union{Nothing, WindFarm}`: Wind farm object containing the initialized simulation state, or `Nothing` if no state is available. See [`WindFarm`](@ref)
- `sim::Sim`: Simulation configuration object containing initialization settings and file paths. See [`Sim`](@ref)

# Returns
- `wf::Union{Nothing, WindFarm}`: The wind farm state, either the original input state (for "init" mode) or a loaded state from disk (for "load" mode)

# Behavior
The function operates in two modes based on `sim.init`:

## "init" Mode
- Uses the provided wind farm state as-is
- If `sim.save_init_state` is `true`, saves the current state to `"T_init.jld2"` in the specified data directory
- Logs the save operation for user feedback

## "load" Mode  
- Attempts to load a previously saved wind farm state from `"T_init.jld2"`
- Falls back to the provided state if loading fails (with warning)
- Handles file I/O errors gracefully

# File Operations
- **Save path**: `\$(sim.path_to_data)/T_init.jld2`
- **Format**: JLD2 binary format for efficient Julia object serialization
- **Error handling**: Loading failures produce warnings but do not halt execution

# Notes
- The function is case-insensitive for the initialization mode string
- File operations use the path specified in `sim.path_to_data`
- Loading errors are caught and logged as warnings, allowing simulation to proceed with the original state
- This mechanism enables reproducible simulations by preserving and reusing initial conditions
"""
function initSimulation(wf::Union{Nothing, WindFarm}, sim::Sim)
    # Initialize the simulation or load an initialized state
    sim_init = lowercase(sim.init)

    if sim_init == "init"
        if sim.save_init_state
            # Save the initialization state to a file
            path = joinpath(sim.path_to_data, "T_init.jld2")
            @info "Saving windfield as $(path) ..."
            jldsave(path; wf)
        end
    elseif sim_init == "load"
        try
            data = load(joinpath(sim.path_to_data, "T_init.jld2"))
            wf = data["wf"]
        catch e
            @warn "Could not load T_init.jld2 from $(sim.path_to_data)\nWill proceed with initialized data." exception=e
        end
    end
    return wf
end

"""
    perturbationOfTheWF!(wf, wind)

Apply stochastic perturbations to the wind field states in-place.

This function adds Gaussian noise to the wind field parameters to model measurement 
uncertainty or natural variability in wind conditions. The perturbations are applied 
conditionally based on the wind perturbation configuration and are added directly 
to the wind farm state matrix.

# Arguments
- `wf`: Wind farm object containing the state matrix `States_WF` to be perturbed
- `wind`: Wind configuration object containing perturbation settings. See [`Wind`](@ref)

# Returns
- `nothing`: The function modifies the wind farm state in-place

# Behavior
The function applies independent Gaussian perturbations to three wind field parameters:

## Velocity Perturbation
- **Condition**: `wind.pertubation.vel == true`
- **Target**: Column 1 of `wf.States_WF` (wind velocity [m/s])
- **Noise**: `wind.pertubation.vel_sigma * randn(nOP × nT)`

## Direction Perturbation  
- **Condition**: `wind.pertubation.dir == true`
- **Target**: Column 2 of `wf.States_WF` (wind direction [degrees])
- **Noise**: `wind.pertubation.dir_sigma * randn(nOP × nT)`

## Turbulence Intensity Perturbation
- **Condition**: `wind.pertubation.ti == true`  
- **Target**: Column 3 of `wf.States_WF` (turbulence intensity [-])
- **Noise**: `wind.pertubation.ti_sigma * randn(nOP × nT)`

# Mathematical Description
For each enabled perturbation type, the function applies:
```
States_WF[:, col] += σ × N(0,1)
```
where:
- `σ` is the standard deviation for the specific parameter
- `N(0,1)` is standard normal random noise with dimensions `(nOP × nT)`
- `nOP` is the number of operational points per turbine
- `nT` is the total number of turbines

# Notes
- The function uses in-place modification (indicated by the `!` suffix)
- Perturbations are applied independently to each operational point and turbine
- The random noise follows a standard normal distribution scaled by the respective sigma values
- Only enabled perturbation types (based on boolean flags) are applied
"""
@views function perturbationOfTheWF!(wf, Wind)
    # perturbationOfTheWF! adds noise to the entire wind field state
    
    # Velocity
    if Wind.pertubation.vel
       wf.States_WF[:, 1] .+= Wind.pertubation.vel_sigma * randn(wf.nOP *wf.nT)
    end

    # Direction
    if Wind.pertubation.dir
       wf.States_WF[:, 2] .+= Wind.pertubation.dir_sigma * randn(wf.nOP *wf.nT)
    end

    # Turbulence Intensity
    if Wind.pertubation.ti
       wf.States_WF[:, 3] .+= Wind.pertubation.ti_sigma * randn(wf.nOP *wf.nT)
    end

    return nothing
end

"""
    findTurbineGroups(wf::WindFarm, floridyn::FloriDyn)

Determine wake interaction dependencies between turbines in a wind farm.

This function analyzes the spatial relationships between turbines to identify which turbines 
are affected by the wakes of upstream turbines. It uses coordinate transformations to the 
wind-aligned reference frame and geometric criteria to determine wake interactions.

# Arguments
- `wf::WindFarm`: Wind farm object containing turbine positions, operational points, and wind field states. See [`WindFarm`](@ref)
- `floridyn::FloriDyn`: FLORIDyn model parameters containing wake interaction thresholds. See [`FloriDyn`](@ref)

# Returns
- `vv_dep::Vector{Vector{Int64}}`: A vector of vectors where `vv_dep[i]` contains the indices of all turbines 
  that affect turbine `i` through wake interactions. Each inner vector lists the upstream turbine indices 
  that influence the wake conditions at the corresponding turbine.

# Algorithm
1. **Coordinate Transformation**: For each turbine pair, transforms coordinates to a wind-aligned frame using the closest operational point
2. **Wake Zone Detection**: Applies geometric criteria to determine if a downstream turbine lies within the wake zone:
   - Upstream extent: `r₁[1] ≥ -uw × D[iT]` (allowing for slight upstream influence)
   - Downstream extent: `r₁[1] ≤ dw × D[iT]` (wake extends downstream)  
   - Lateral extent: `|r₁[2]| ≤ cw × D[iT]` (wake width constraint)
3. **Dependency Matrix**: Constructs a boolean dependency matrix and extracts indices for each turbine

# Mathematical Description
The wake interaction criteria are evaluated in the wind-aligned coordinate system:
```
r₁ = R(φ) × (rₒₚ - rₜᵤᵣᵦ)
```
where:
- `R(φ)` is the rotation matrix for wind direction angle `φ`
- `rₒₚ` is the position of the closest operational point from the upstream turbine
- `rₜᵤᵣᵦ` is the position of the downstream turbine being evaluated

# Notes
- The function uses the closest operational point from each upstream turbine to determine wind direction
- Wake zones are defined as multiples of rotor diameter using the FLORIDyn parameters
- Self-interaction (turbine affecting itself) is explicitly excluded
- The coordinate transformation accounts for the SOWFA wind direction convention
"""
@views function findTurbineGroups(wf::WindFarm, floridyn::FloriDyn)
    # Extract parameters from settings struct
    dw = floridyn.deltaDW
    cw = floridyn.deltaCW
    uw = floridyn.deltaUW

    # Initialize outputs
    vv_dep = Vector{Vector{Int64}}(undef,wf.nT)  # Equivalent of cell array in MATLAB[1][3]
    dep  = falses(wf.nT, wf.nT)

    for iT in 1:wf.nT
        for iiT in 1:wf.nT
            if iiT == iT
                continue
            end

            # Find closest OP from wake of iT to turbine iiT without allocations
            start_idx = wf.StartI[iT]
            end_idx = start_idx + wf.nOP - 1
            
            min_dist_sq = Inf
            I_op = start_idx
            
            # Manual loop to find minimum distance without allocations
            @inbounds for op_idx in start_idx:end_idx
                dx = wf.States_OP[op_idx, 1] - wf.posBase[iiT, 1]
                dy = wf.States_OP[op_idx, 2] - wf.posBase[iiT, 2]
                dist_sq = dx * dx + dy * dy
                
                if dist_sq < min_dist_sq
                    min_dist_sq = dist_sq
                    I_op = op_idx
                end
            end

            # Compute angle and relative vector without allocations
            phi = angSOWFA2world(wf.States_WF[I_op, 2])
            cos_phi = cos(phi)
            sin_phi = sin(phi)
            
            # Relative vector r0 = States_OP - posBase
            r0_x = wf.States_OP[I_op, 1] - wf.posBase[iiT, 1]
            r0_y = wf.States_OP[I_op, 2] - wf.posBase[iiT, 2]
            
            # Apply rotation matrix R01(phi) * r0 manually
            r1_x = cos_phi * r0_x + sin_phi * r0_y
            r1_y = -sin_phi * r0_x + cos_phi * r0_y

            # Apply dependency check
            if (-r1_x <= uw*wf.D[iT]) && (r1_x <= dw*wf.D[iT]) && (abs(r1_y) <= cw*wf.D[iT])
                dep[iiT, iT] = true
            end
        end
    end

    # Extract indices manually to avoid allocations in findall
    for iT in 1:wf.nT
        # Count true values first
        count = 0
        @inbounds for iiT in 1:wf.nT
            if dep[iT, iiT]
                count += 1
            end
        end
        
        # Allocate vector of exact size and fill it
        vv_dep[iT] = Vector{Int64}(undef, count)
        idx = 1
        @inbounds for iiT in 1:wf.nT
            if dep[iT, iiT]
                vv_dep[iT][idx] = iiT
                idx += 1
            end
        end
    end

    return vv_dep
end

"""
    interpolateOPs(wf)

Compute interpolation weights and indices for operational points affecting each turbine.

This function determines the optimal interpolation strategy for each turbine by identifying 
the closest operational points from influencing upstream turbines. It computes weights and 
indices that enable smooth interpolation of wind field states and turbine conditions at 
arbitrary turbine positions.

# Arguments
- `wf`: Wind farm object containing turbine dependencies, operational point states, and positional data
  - `wf.nT`: Number of turbines
  - `wf.StartI`: Starting indices for each turbine's operational points  
  - `wf.dep`: Dependency relationships between turbines (from [`findTurbineGroups`](@ref))
  - `wf.States_OP`: Matrix of operational point states
  - `wf.posBase`: Base positions of turbines [m]
  - `wf.nOP`: Number of operational points per turbine

# Returns
- `intOPs::Vector{Matrix{Float64}}`: Interpolation data for each turbine where `intOPs[i]` is an 
  `N×4` matrix for turbine `i` with `N` influencing turbines. Each row contains:
  - Column 1: First operational point index
  - Column 2: Weight for first operational point
  - Column 3: Second operational point index  
  - Column 4: Weight for second operational point

# Algorithm
For each turbine and its influencing upstream turbines:

1. **Distance Calculation**: Computes Euclidean distances from all operational points of 
   the influencing turbine to the target turbine position

2. **Interpolation Strategy Selection**: Based on the closest operational point location:
   - **First OP closest**: Uses first and second operational points
   - **Last OP closest**: Uses second-to-last and last operational points
   - **Interior OP closest**: Uses the two closest operational points for optimal interpolation

3. **Weight Computation**: For interior cases, applies linear projection to determine interpolation weights:
   ```julia
   d = dot(ab, ac) / dot(ab, ab)
   weights = [1-d, d] # Clamped to [0,1]
   ```

# Mathematical Description
The interpolation uses linear projection for weight computation:
```
d = (b - a) · (c - a) / |b - a|²
```
where:
- `a`, `b` are positions of the two closest operational points
- `c` is the target turbine position
- `d` is the projection parameter (clamped to [0,1])

# Notes
- Edge cases (first/last operational points) use predefined weight combinations
- Weights always sum to 1.0 for proper interpolation
- The function handles variable numbers of influencing turbines per target turbine
- Interpolation indices are global across the entire operational point matrix
- This preprocessing enables efficient interpolation during simulation time steps
"""
function interpolateOPs(wf)
    @assert length(wf.dep) > 0 "No dependencies found! Ensure `findTurbineGroups` was called first."
    intOPs = Vector{Matrix{Float64}}(undef,wf.nT)  # Cell equivalent in Julia

    for iT in 1:wf.nT  # For every turbine
        intOPs[iT] = zeros(length(wf.dep[iT]), 4)

        for iiT in 1:length(wf.dep[iT])  # for every influencing turbine
            iiaT = wf.dep[iT][iiT]       # actual turbine index

            # Compute distances from OPs of turbine iiaT to current turbine
            start_idx    = wf.StartI[iiaT]
            OP_positions = wf.States_OP[start_idx:(start_idx +wf.nOP - 1), 1:2]
            turb_pos     = wf.posBase[iT, 1:2]

            # Euclidean distances to the turbine position
            dist = sqrt.(sum((OP_positions' .- turb_pos).^2, dims=2))
            dist = vec(dist)  # make it a flat vector

            # Indices of sorted distances
            sorted_indices = sortperm(dist)

            if sorted_indices[1] == 1
                # Closest is first OP (unlikely)
                intOPs[iT][iiT, :] = [wf.StartI[iiaT], 1.0,wf.StartI[iiaT] + 1, 0.0]
            elseif sorted_indices[1] == wf.nOP
                # Closest is last OP (possible)
                intOPs[iT][iiT, :] = [wf.StartI[iiaT] + wf.nOP - 2, 0.0, wf.StartI[iiaT] + wf.nOP - 1, 1.0]
            else
                # Use two closest OPs for interpolation
                indOP1 = wf.StartI[iiaT] - 1 + sorted_indices[1]
                indOP2 = wf.StartI[iiaT] - 1 + sorted_indices[2]

                a = wf.States_OP[indOP1, 1:2]
                b = wf.States_OP[indOP2, 1:2]
                c = wf.posBase[iT, 1:2]

                ab = b .- a
                ac = c .- a
                d  = dot(ab, ac) / dot(ab, ab)
                d  = clamp(d, 0.0, 1.0)

                r1 = 1.0 - d
                r2 = d

                intOPs[iT][iiT, :] = [indOP1, r1, indOP2, r2]
            end
        end
    end

    return intOPs
end

"""
    setUpTmpWFAndRun(set::Settings, wf, floris::Floris, wind::Wind)

Execute FLORIS wake calculations for all turbines in a wind farm with wake interactions.

This function orchestrates the computation of wake effects for each turbine by setting up 
temporary wind farm configurations that include influencing upstream turbines. It handles 
both single turbine (no wake interactions) and multi-turbine scenarios with complex wake 
interaction patterns.

# Arguments
- `set::Settings`: Simulation settings and configuration parameters
- `wf`: Wind farm object containing turbine positions, operational points, dependencies, and interpolation data
  - `wf.nT`: Number of turbines
  - `wf.States_WF`: Wind field states matrix
  - `wf.States_T`: Turbine states matrix
  - `wf.States_OP`: Operational point states matrix
  - `wf.dep`: Turbine dependency relationships (from [`findTurbineGroups`](@ref))
  - `wf.intOPs`: Interpolation weights and indices (from [`interpolateOPs`](@ref))
  - `wf.posBase`: Base turbine positions [m]
  - `wf.posNac`: Nacelle position offsets [m]
  - `wf.D`: Rotor diameters [m]
  - `wf.StartI`: Starting indices for each turbine's operational points
- `floris::Floris`: FLORIS model parameters for wake calculations. See [`Floris`](@ref)
- `wind::Wind`: Wind field configuration including shear properties. See [`Wind`](@ref)

# Returns
- `M::Matrix{Float64}`: Results matrix of size `(nT × 3)` where each row contains:
  - Column 1: Total velocity reduction factor (product of all wake effects)
  - Column 2: Combined added turbulence intensity from all upstream turbines
  - Column 3: Effective wind speed at turbine [m/s]
- `wf`: Updated wind farm object with modified fields:
  - `wf.Weight`: Normalized interpolation weights for each turbine
  - `wf.red_arr`: Wake reduction matrix showing turbine-to-turbine wake effects

# Algorithm
The function processes each turbine individually:

## Single Turbine Case (No Dependencies)
- Directly calls FLORIS with the turbine's wind field state
- No wake interactions considered
- Results stored directly in output matrix

## Multi-Turbine Case (With Dependencies)  
1. **Temporary Configuration Setup**: Creates temporary arrays sized for the target turbine plus all influencing turbines
2. **Interpolation Application**: Uses precomputed interpolation weights to determine states at influencing turbine positions
3. **Coordinate Transformation**: Applies wind direction-based coordinate transformations to account for spatial offsets
4. **FLORIS Execution**: Runs wake model with the complete multi-turbine configuration
5. **Result Processing**: Combines wake effects and normalizes weights

# Mathematical Description
For multi-turbine scenarios, the effective position of influencing turbines is computed as:
```
tmp_Tpos[i] = base_position - R(φ) × [offset_x, offset_y, offset_z]
```
where `R(φ)` is the rotation matrix for wind direction `φ`.

The total wake reduction is the product of individual wake effects:
```
T_red = ∏ᵢ T_red_arr[i]
```

Combined turbulence intensity follows root-sum-square combination:
```
T_addedTI = √(∑ᵢ T_aTI_arr[i]²)
```

# Wind Field Interpolation
The function supports optional wind field interpolation via coefficient matrices:
- **Velocity interpolation**: Uses `wf.C_Vel` if available
- **Direction interpolation**: Uses `wf.C_Dir` if available

# Notes
- The function modifies the wind farm object in-place, updating weight and reduction arrays
- Interpolation weights are normalized to ensure proper weighting
- Special handling for variable rotor diameter configurations
- Coordinate transformations use the SOWFA to world conversion via [`angSOWFA2world`](@ref)
- The algorithm efficiently handles both simple single-turbine and complex multi-turbine wake scenarios
"""
@views function setUpTmpWFAndRun(set::Settings, wf, floris::Floris, wind::Wind)
    # Initialize outputs
    M = zeros(wf.nT, 3)
    wf.Weight = Vector{Vector{Float64}}(undef,wf.nT)
    wf.red_arr = ones(wf.nT,wf.nT)

    for iT in 1:wf.nT
        # Interpolate Wind field if needed
        iTWFState = copy(wf.States_WF[wf.StartI[iT], :])

        if hasfield(typeof(wf), :C_Vel)
            iTWFState[1] = dot(wf.C_Vel[iT, :],wf.States_WF[:, 1])
        end

        if hasfield(typeof(wf), :C_Dir)
            iTWFState[2] = dot(wf.C_Dir[iT, :],wf.States_WF[:, 2])
        end

        if isempty(wf.dep[iT])
            # Single turbine case
            T_red_arr, _, _ = runFLORIS(
                set,
                (wf.posBase[iT,:] +wf.posNac[iT,:])',
                iTWFState',
               wf.States_T[wf.StartI[iT], :]',
               wf.D[iT],
                floris,
                wind.shear
            )
            M[iT, :] = [T_red_arr, 0, T_red_arr *wf.States_WF[wf.StartI[iT], 1]]
           wf.red_arr[iT, iT] = T_red_arr
            continue
        end

        # Multi-turbine setup
        tmp_nT = length(wf.dep[iT]) + 1

        tmp_Tpos = repeat(wf.posBase[iT,:]' + wf.posNac[iT,:]', tmp_nT)
        tmp_WF   = repeat(iTWFState', tmp_nT)
        tmp_Tst  = repeat((wf.States_T[wf.StartI[iT], :])', tmp_nT)

        tmp_D = if wf.D[end] > 0
            vcat(wf.D[wf.dep[iT]],wf.D[iT])
        else
           wf.D
        end

        for iiT in 1:(tmp_nT - 1)
            OP1_i = Int(wf.intOPs[iT][iiT, 1])  # Index OP 1
            OP1_r = wf.intOPs[iT][iiT, 2]       # Ratio OP 1
            OP2_i = Int(wf.intOPs[iT][iiT, 3])  # Index OP 2
            OP2_r = wf.intOPs[iT][iiT, 4]       # Ratio OP 2

            OPi_l = OP1_r * wf.States_OP[OP1_i, :] + OP2_r * wf.States_OP[OP2_i, :]
            tmp_Tpos[iiT, :] = OPi_l[1:3]
            tmp_Tst[iiT, :] = OP1_r *wf.States_T[OP1_i, :] + OP2_r *wf.States_T[OP2_i, :]
            tmp_WF[iiT, :]  = OP1_r *wf.States_WF[OP1_i, :] + OP2_r *wf.States_WF[OP2_i, :]

            si = wf.StartI[wf.dep[iT][iiT]]

            if hasfield(typeof(wf), :C_Vel)
                C_weights = wf.C_Vel[iT, si:(si + wf.nOP - 1)]
                C_weights ./= sum(C_weights)
                tmp_WF[iiT, 1] = dot(C_weights, wf.States_WF[si:si + wf.nOP - 1, 1])
            end
            if hasfield(typeof(wf), :C_Dir)
                C_weights = wf.C_Dir[iT, si:(si + wf.nOP - 1)]
                C_weights ./= sum(C_weights)
                tmp_WF[iiT, 2] = dot(C_weights, wf.States_WF[si:si + wf.nOP - 1, 2])
            end

            tmp_phi = size(tmp_WF, 2) == 4 ? angSOWFA2world(tmp_WF[iiT, 4]) : angSOWFA2world(tmp_WF[iiT, 2])

            tmp_Tpos[iiT, 1] -= cos(tmp_phi) * OPi_l[4] - sin(tmp_phi) * OPi_l[5]
            tmp_Tpos[iiT, 2] -= sin(tmp_phi) * OPi_l[4] + cos(tmp_phi) * OPi_l[5]
            tmp_Tpos[iiT, 3] -= OPi_l[6]
        end

        # Run FLORIS                
        T_red_arr, T_aTI_arr, T_Ueff, T_weight = runFLORIS(set, tmp_Tpos, tmp_WF, tmp_Tst, tmp_D, floris, wind.shear)

        T_red = prod(T_red_arr)
        wf.red_arr[iT, vcat(wf.dep[iT], iT)] = T_red_arr
        T_addedTI = sqrt(sum(T_aTI_arr .^ 2))
        wf.Weight[iT] = T_weight

        if wf.D[end] <= 0
            dists = zeros(tmp_nT - 1)
            plot_WF = zeros(tmp_nT - 1, size(wf.States_WF, 2))
            plot_OP = zeros(tmp_nT - 1, 2)
            for iiT in 1:(tmp_nT - 1)
                OP1_i, OP1_r, OP2_i, OP2_r =wf.intOPs[iT][iiT, :]
                OPi_l = OP1_r *wf.States_OP[OP1_i, :] + OP2_r *wf.States_OP[OP2_i, :]
                plot_OP[iiT, :] = OPi_l[1:2]
                plot_WF[iiT, :] = OP1_r *wf.States_WF[OP1_i, :] + OP2_r *wf.States_WF[OP2_i, :]
                dists[iiT] = norm(OPi_l[1:2] .-wf.posBase[iT,1:2])
            end

            I = sortperm(dists)
            if length(I) == 1
                Ufree = plot_WF[I[1], 1]
                T_Ueff = T_red * Ufree
            else
                a = plot_OP[I[1], :]'
                b = plot_OP[I[2], :]'
                c =wf.posBase[iT, 1:2]'
                d = clamp((dot(b - a, c - a)) / dot(b - a, b - a), 0.0, 1.0)
                r1, r2 = 1.0 - d, d
                Ufree = r1 * plot_WF[I[1], 1] + r2 * plot_WF[I[2], 1]
                T_Ueff = T_red * Ufree
            end
        end

        M[iT, :] = [T_red, T_addedTI, T_Ueff]

        wS = sum(wf.Weight[iT])
        if wS > 0
           wf.Weight[iT] =wf.Weight[iT] ./ wS
        else
           wf.Weight[iT] .= 0.0
        end
    end

    return M, wf
end


"""
    runFLORIDyn(set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, 
                floridyn::FloriDyn, floris::Floris)

Main entry point for the FLORIDyn closed-loop simulation.

# Arguments
- `set::Settings`: Simulation settings and configuration parameters.
- `wf::WindFarm`: See: [WindFarm](@ref) simulation state, including turbine and wind farm states.
- `wind::Wind`: See: [Wind](@ref) field settings.
- `sim::Sim`: Simulation state or configuration object. See: [`Sim`](@ref)
- `con::Con`: Controller object or control parameters. See: [`Con`](@ref)
- `floridyn::FloriDyn`: Parameters specific to the FLORIDyn model. See: [`FloriDyn`](@ref)
- `floris::Floris`: Parameters specific to the FLORIS model. See: [`Floris`](@ref)

# Returns
A tuple `(wf, md, mi)` containing:
- `wf::WindFarm`: Updated simulation state with final turbine positions, wind field states, and operational point data
- `md::DataFrame`: Measurement data with columns:
  - `:Time`: Simulation time steps
  - `:ForeignReduction`: Wind speed reduction factors (%) due to wake effects from other turbines
  - `:AddedTurbulence`: Additional turbulence intensity (%) induced by upstream turbines
  - `:EffWindSpeed`: Effective wind speed (m/s) at each turbine after wake effects
  - `:FreeWindSpeed`: Free-stream wind speed (m/s) without wake interference
  - `:PowerGen`: Generated electrical power (MW) for each turbine
- `mi::Matrix`: Interaction matrix combining time data with turbine-to-turbine wake interaction coefficients 
                for each simulation step

# Description
Runs a closed-loop wind farm simulation using the FLORIDyn and FLORIS models, 
applying control strategies and updating turbine states over time.

"""
function runFLORIDyn(set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, floridyn::FloriDyn, floris::Floris)
    nT      = wf.nT
    nSim    = sim.n_sim_steps
    M       = zeros(nSim * nT, 6)
    M[:, 1] .= 1.0  # Set first column to 1
    M_int   = Vector{Matrix{Float64}}(undef, nSim)

    SimTime = sim.start_time
    
    buffers = FLORIDyn.IterateOPsBuffers(wf)
    for it in 1:nSim
        sim.sim_step = it

        # ========== PREDICTION ==========
        iterateOPs!(set.iterate_mode, wf, sim, floris, floridyn, buffers)

        # ========== Wind Field Perturbation ==========
        perturbationOfTheWF!(wf, wind)

        # ========== Get FLORIS reductions ==========
        wf.dep = findTurbineGroups(wf, floridyn)
        wf.intOPs = interpolateOPs(wf)
        a, b = setUpTmpWFAndRun(set, wf, floris, wind)
        tmpM, wf = a, b
        M[(it-1)*nT+1 : it*nT, 2:4] .= tmpM
        M[(it-1)*nT+1 : it*nT, 1]   .= SimTime
        wf.States_T[wf.StartI, 3] = tmpM[:, 2]
        M_int[it] = wf.red_arr

        # ========== wind field corrections ==========
        wf, wind = correctVel(set.cor_vel_mode, set, wf, wind, SimTime, floris, tmpM)
        correctDir!(set.cor_dir_mode, set, wf, wind, SimTime)
        correctTI!(set.cor_turb_mode, set, wf, wind, SimTime)

        # Save free wind speed as measurement
        M[(it-1)*nT+1 : it*nT, 5] = wf.States_WF[wf.StartI, 1]

        # ========== Get Control settings ==========
        wf.States_T[wf.StartI, 2] = (
            wf.States_WF[wf.StartI, 2] .-
                getYaw(set.control_mode, con.yaw_data, (1:nT), SimTime)'
        )

        # ========== Calculate Power ==========
        P = getPower(wf, tmpM, floris, con)
        M[(it-1)*nT+1:it*nT, 6] = P

        SimTime += sim.time_step
    end
    # Convert `M` to DataFrame and scale measurements
    md = DataFrame(
        (M * diagm([1; 100; 100; 1; 1; 1e-6])),
        [:Time, :ForeignReduction, :AddedTurbulence, :EffWindSpeed, :FreeWindSpeed, :PowerGen]
    )
    mi = hcat(md.Time, hcat(M_int...)')
    return wf, md, mi
end
