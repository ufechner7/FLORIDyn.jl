# Copyright (c) 2025 Marcus Becker, Uwe Fechner
# SPDX-License-Identifier: BSD-3-Clause

"""
    angSOWFA2world(deg_SOWFA) -> Float64

Convert wind direction angle from SOWFA convention to world coordinate system.

This function performs coordinate transformation between different wind direction conventions
used in wind farm simulations. SOWFA (Simulator fOr Wind Farm Applications) uses a different
angular reference system than the standard world coordinate system used in calculations.

# Arguments
- `deg_SOWFA::Real`: Wind direction angle in SOWFA convention [degrees]

# Returns
- `rad_World`: Wind direction angle in world coordinate system [radians]

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
    initSimulation(wf::Union{Nothing, WindFarm}, sim::Sim) -> Union{Nothing, WindFarm}

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
    perturbationOfTheWF!(wf::WindFarm, wind::Wind) -> Nothing

Apply stochastic perturbations to the wind field states in-place.

This function adds Gaussian noise to the wind field parameters to model measurement 
uncertainty or natural variability in wind conditions. The perturbations are applied 
conditionally based on the wind perturbation configuration and are added directly 
to the wind farm state matrix.

# Arguments
- `wf::WindFarm`: Wind farm struct containing the state matrix `States_WF` to be perturbed
- `wind::Wind`: Wind configuration struct containing perturbation settings. See [`Wind`](@ref)

# Returns
- `nothing`: The function modifies the wind farm state in-place

# Behavior
The function applies independent Gaussian perturbations to three wind field parameters:

## Velocity Perturbation
- **Condition**: `wind.perturbation.vel == true`
- **Target**: Column 1 of `wf.States_WF` (wind velocity [m/s])
- **Noise**: `wind.perturbation.vel_sigma * randn(nOP × nT)`

## Direction Perturbation  
- **Condition**: `wind.perturbation.dir == true`
- **Target**: Column 2 of `wf.States_WF` (wind direction [degrees])
- **Noise**: `wind.perturbation.dir_sigma * randn(nOP × nT)`

## Turbulence Intensity Perturbation
- **Condition**: `wind.perturbation.ti == true`  
- **Target**: Column 3 of `wf.States_WF` (turbulence intensity [-])
- **Noise**: `wind.perturbation.ti_sigma * randn(nOP × nT)`

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
@views function perturbationOfTheWF!(wf::WindFarm, wind::Wind)
    # perturbationOfTheWF! adds noise to the entire wind field state
    
    # Velocity
    if wind.perturbation.vel
       wf.States_WF[:, 1] .+= wind.perturbation.vel_sigma * randn(wf.nOP *wf.nT)
    end

    # Direction
    if wind.perturbation.dir
       wf.States_WF[:, 2] .+= wind.perturbation.dir_sigma * randn(wf.nOP *wf.nT)
    end

    # Turbulence Intensity
    if wind.perturbation.ti
       wf.States_WF[:, 3] .+= wind.perturbation.ti_sigma * randn(wf.nOP *wf.nT)
    end

    return nothing
end

"""
    findTurbineGroups(wf::WindFarm, floridyn::FloriDyn) -> Vector{Vector{Int64}}

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
    interpolateOPs!(unified_buffers::UnifiedBuffers, intOPs::Vector{Matrix{Float64}}, wf::WindFarm)

Compute interpolation weights and indices for operational points affecting each turbine using a unified buffer.

This function performs interpolation calculations while avoiding 
memory allocations by reusing pre-allocated buffer arrays from a unified buffer struct. 
This is critical for performance when called repeatedly in loops, such as in flow field calculations.

# Arguments
- `intOPs::Vector{Matrix{Float64}}`: Pre-allocated vector of matrices to store interpolation results
- `wf::WindFarm`: Wind farm object containing turbine positions and operational point data
- `unified_buffers::UnifiedBuffers`: Unified buffer struct containing pre-allocated arrays including:
  - `dist_buffer`: Buffer for distance calculations (length ≥ wf.nOP)
  - `sorted_indices_buffer`: Buffer for sorting indices (length ≥ wf.nOP)

# Returns
- `intOPs::Vector{Matrix{Float64}}`: Results filled in-place

# Performance Notes
- All temporary arrays are reused from the unified buffer struct
- No memory allocations occur during execution
- Suitable for use in hot loops and parallel contexts

# Example
```julia
# Create unified buffers
unified_buffers = create_unified_buffers(wf)

# Pre-allocate interpolation matrices
intOPs = [zeros(length(wf.dep[iT]), 4) for iT in 1:wf.nT]

# Non-allocating interpolation
interpolateOPs!(unified_buffers, intOPs, wf)
```
"""
function interpolateOPs!(unified_buffers::UnifiedBuffers, intOPs::Vector{Matrix{Float64}}, wf::WindFarm)
    @assert length(wf.dep) > 0 "No dependencies found! Ensure `findTurbineGroups` was called first."

    for iT in 1:wf.nT  # For every turbine
        for iiT in 1:length(wf.dep[iT])  # for every influencing turbine
            iiaT = wf.dep[iT][iiT]       # actual turbine index

            # Compute distances from OPs of turbine iiaT to current turbine
            start_idx = wf.StartI[iiaT]
            turb_pos_x = wf.posBase[iT, 1]
            turb_pos_y = wf.posBase[iT, 2]

            # Compute Euclidean distances to the turbine position (non-allocating)
            @views dist = unified_buffers.dist_buffer[1:wf.nOP]
            for i in 1:wf.nOP
                op_idx = start_idx + i - 1
                dx = wf.States_OP[op_idx, 1] - turb_pos_x
                dy = wf.States_OP[op_idx, 2] - turb_pos_y
                dist[i] = sqrt(dx * dx + dy * dy)
            end

            # Sort indices by distance (reuse buffer)
            @views sorted_indices = unified_buffers.sorted_indices_buffer[1:wf.nOP]
            for i in 1:wf.nOP
                sorted_indices[i] = i
            end
            sort!(sorted_indices, by=i -> dist[i])

            if sorted_indices[1] == 1
                # Closest is first OP (unlikely)
                intOPs[iT][iiT, 1] = wf.StartI[iiaT]
                intOPs[iT][iiT, 2] = 1.0
                intOPs[iT][iiT, 3] = wf.StartI[iiaT] + 1
                intOPs[iT][iiT, 4] = 0.0
            elseif sorted_indices[1] == wf.nOP
                # Closest is last OP (possible)
                intOPs[iT][iiT, 1] = wf.StartI[iiaT] + wf.nOP - 2
                intOPs[iT][iiT, 2] = 0.0
                intOPs[iT][iiT, 3] = wf.StartI[iiaT] + wf.nOP - 1
                intOPs[iT][iiT, 4] = 1.0
            else
                # Use two closest OPs for interpolation
                indOP1 = wf.StartI[iiaT] - 1 + sorted_indices[1]
                indOP2 = wf.StartI[iiaT] - 1 + sorted_indices[2]

                a_x = wf.States_OP[indOP1, 1]
                a_y = wf.States_OP[indOP1, 2]
                b_x = wf.States_OP[indOP2, 1]
                b_y = wf.States_OP[indOP2, 2]
                c_x = wf.posBase[iT, 1]
                c_y = wf.posBase[iT, 2]

                ab_x = b_x - a_x
                ab_y = b_y - a_y
                ac_x = c_x - a_x
                ac_y = c_y - a_y
                
                d = (ab_x * ac_x + ab_y * ac_y) / (ab_x * ab_x + ab_y * ab_y)
                d = clamp(d, 0.0, 1.0)

                r1 = 1.0 - d
                r2 = d

                intOPs[iT][iiT, 1] = indOP1
                intOPs[iT][iiT, 2] = r1
                intOPs[iT][iiT, 3] = indOP2
                intOPs[iT][iiT, 4] = r2
            end
        end
    end

    return intOPs
end

"""
    setUpTmpWFAndRun(set::Settings, wf::WindFarm, 
                     floris::Floris, wind::Wind) --> (Matrix, WindFarm)

Execute FLORIS wake calculations for all turbines in a wind farm with wake interactions.

This function orchestrates the computation of wake effects for each turbine by setting up 
temporary wind farm configurations that include influencing upstream turbines. It handles 
both single turbine (no wake interactions) and multi-turbine scenarios with complex wake 
interaction patterns.

# Arguments
- `set::Settings`: Simulation settings and configuration parameters
- `wf::WindFarm`: Wind farm object containing turbine positions, operational points, dependencies, and interpolation data
  - `wf.nT`: Number of turbines
  - `wf.States_WF`: Wind field states matrix
  - `wf.States_T`: Turbine states matrix
  - `wf.States_OP`: Operational point states matrix
  - `wf.dep`: Turbine dependency relationships (from [`findTurbineGroups`](@ref))
  - `wf.intOPs`: Interpolation weights and indices (from [`interpolateOPs!`](@ref))
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
- `wf::WindFarm`: Updated wind farm object with modified fields:
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
@views function setUpTmpWFAndRun(set::Settings, wf::WindFarm, floris::Floris, wind::Wind)
    # Initialize outputs
    M = zeros(wf.nT, 3)
    wf.Weight = [Float64[] for _ in 1:wf.nT]
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
                OP1_i_f, OP1_r, OP2_i_f, OP2_r = wf.intOPs[iT][iiT, :]
                OP1_i = Int(round(OP1_i_f))
                OP2_i = Int(round(OP2_i_f))
                OPi_l = OP1_r * wf.States_OP[OP1_i, :] + OP2_r * wf.States_OP[OP2_i, :]
                plot_OP[iiT, :] = OPi_l[1:2]
                plot_WF[iiT, :] = OP1_r * wf.States_WF[OP1_i, :] + OP2_r * wf.States_WF[OP2_i, :]
                dists[iiT] = norm(OPi_l[1:2] .- wf.posBase[iT,1:2])
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
    setUpTmpWFAndRun!(ub::UnifiedBuffers, wf::WindFarm, set::Settings, floris::Floris, wind::Wind) 
                      -> (Matrix{Float64}, WindFarm)

Non-allocating version of [`setUpTmpWFAndRun`](@ref) that uses a unified buffer struct.

This function performs the same calculations as `setUpTmpWFAndRun` but avoids memory allocations
by reusing pre-allocated buffer arrays from a [`UnifiedBuffers`](@ref) struct. This is particularly 
important for parallel execution and performance-critical loops where garbage collection overhead 
needs to be minimized.

# Arguments
- `ub::UnifiedBuffers`: Unified buffer struct containing all pre-allocated arrays
  - `ub.M_buffer`: Pre-allocated buffer for results matrix (size: nT × 3)
  - `ub.iTWFState_buffer`: Buffer for turbine wind field state
  - `ub.tmp_Tpos_buffer`: Buffer for temporary turbine positions
  - `ub.tmp_WF_buffer`: Buffer for temporary wind field states
  - `ub.tmp_Tst_buffer`: Buffer for temporary turbine states
  - `ub.dists_buffer`: Buffer for distance calculations
  - `ub.plot_WF_buffer`: Buffer for plotting wind field data
  - `ub.plot_OP_buffer`: Buffer for plotting operating point data
- `wf::WindFarm`: Wind farm object containing turbine data
- `set::Settings`: Settings object containing simulation parameters
- `floris::Floris`: FLORIS model parameters for wake calculations
- `wind::Wind`: Wind field configuration

# Returns
- `M::Matrix{Float64}`: Same as the input `ub.M_buffer`, filled with results
- `wf::WindFarm`: Modified wind farm object with updated internal state

# Performance Notes
- Uses in-place operations to minimize memory allocations
- Buffers must be pre-sized correctly for the specific wind farm configuration
- Thread-safe when each thread uses its own set of buffers
"""
function setUpTmpWFAndRun!(ub::UnifiedBuffers, wf::WindFarm, set::Settings, floris::Floris, wind::Wind; alloc=nothing)
    # Reuse the provided M_buffer instead of allocating new
    ub.M_buffer .= 0.0  # Clear the buffer
    wf.Weight = [Float64[] for _ in 1:wf.nT]
    wf.red_arr = ones(wf.nT, wf.nT)  # Initialize if not already allocated

    if !isnothing(alloc)
        alloc.m += 1
    end


    a = @allocated for iT in 1:wf.nT
        # Reuse iTWFState_buffer instead of allocating
        ub.iTWFState_buffer .= wf.States_WF[wf.StartI[iT], :]

        if hasfield(typeof(wf), :C_Vel)
            ub.iTWFState_buffer[1] = dot(wf.C_Vel[iT, :],wf.States_WF[:, 1])
        end

        if hasfield(typeof(wf), :C_Dir)
            ub.iTWFState_buffer[2] = dot(wf.C_Dir[iT, :],wf.States_WF[:, 2])
        end

        # @assert !isempty(wf.dep[iT])

        if isempty(wf.dep[iT])
            # Single turbine case
            T_red_arr, _, _ = runFLORIS(
                set,
                (wf.posBase[iT,:] +wf.posNac[iT,:])',
                ub.iTWFState_buffer',
               wf.States_T[wf.StartI[iT], :]',
               wf.D[iT],
                floris,
                wind.shear
            )
            ub.M_buffer[iT, :] = [T_red_arr, 0, T_red_arr * wf.States_WF[wf.StartI[iT], 1]]
            wf.red_arr[iT, iT] = T_red_arr
            continue
        end

        # Multi-turbine setup using pre-allocated buffers
        tmp_nT = length(wf.dep[iT]) + 1

        # Reuse buffers instead of repeat operations
        a = @allocated for row in 1:tmp_nT
            ub.tmp_Tpos_buffer[row, :] = wf.posBase[iT,:]' + wf.posNac[iT,:]'
            ub.tmp_WF_buffer[row, :] = ub.iTWFState_buffer'
            ub.tmp_Tst_buffer[row, :] = (wf.States_T[wf.StartI[iT], :])'
        end
        if !isnothing(alloc)
            alloc.for2 += a
        end

        a = @allocated tmp_D = if wf.D[end] > 0
            vcat(wf.D[wf.dep[iT]], wf.D[iT])
        else
           wf.D
        end

        if !isnothing(alloc)
            alloc.if1 += a
        end

        # @assert tmp_nT > 1

        @allocated for iiT in 1:(tmp_nT - 1)
            OP1_i = Int(wf.intOPs[iT][iiT, 1])  # Index OP 1
            OP1_r = wf.intOPs[iT][iiT, 2]       # Ratio OP 1
            OP2_i = Int(wf.intOPs[iT][iiT, 3])  # Index OP 2
            OP2_r = wf.intOPs[iT][iiT, 4]       # Ratio OP 2

            OPi_l = OP1_r * wf.States_OP[OP1_i, :] + OP2_r * wf.States_OP[OP2_i, :]
            ub.tmp_Tpos_buffer[iiT, :] = OPi_l[1:3]
            ub.tmp_Tst_buffer[iiT, :] = OP1_r *wf.States_T[OP1_i, :] + OP2_r *wf.States_T[OP2_i, :]
            ub.tmp_WF_buffer[iiT, :]  = OP1_r *wf.States_WF[OP1_i, :] + OP2_r *wf.States_WF[OP2_i, :]

            si = wf.StartI[wf.dep[iT][iiT]]

            if hasfield(typeof(wf), :C_Vel)
                C_weights = wf.C_Vel[iT, si:(si + wf.nOP - 1)]
                C_weights ./= sum(C_weights)
                ub.tmp_WF_buffer[iiT, 1] = dot(C_weights, wf.States_WF[si:si + wf.nOP - 1, 1])
            end
            if hasfield(typeof(wf), :C_Dir)
                C_weights = wf.C_Dir[iT, si:(si + wf.nOP - 1)]
                C_weights ./= sum(C_weights)
                ub.tmp_WF_buffer[iiT, 2] = dot(C_weights, wf.States_WF[si:si + wf.nOP - 1, 2])
            end

            tmp_phi = size(ub.tmp_WF_buffer, 2) == 4 ? angSOWFA2world(ub.tmp_WF_buffer[iiT, 4]) : angSOWFA2world(ub.tmp_WF_buffer[iiT, 2])

            ub.tmp_Tpos_buffer[iiT, 1] -= cos(tmp_phi) * OPi_l[4] - sin(tmp_phi) * OPi_l[5]
            ub.tmp_Tpos_buffer[iiT, 2] -= sin(tmp_phi) * OPi_l[4] + cos(tmp_phi) * OPi_l[5]
            ub.tmp_Tpos_buffer[iiT, 3] -= OPi_l[6]
        end

        if !isnothing(alloc)
            alloc.for3 += a
        end

        # Run FLORIS using the buffer views
        tmp_Tpos_view = @view ub.tmp_Tpos_buffer[1:tmp_nT, :]
        tmp_WF_view = @view ub.tmp_WF_buffer[1:tmp_nT, :]
        tmp_Tst_view = @view ub.tmp_Tst_buffer[1:tmp_nT, :]
        a = @allocated T_red_arr, T_aTI_arr, T_Ueff, T_weight = runFLORIS(set, tmp_Tpos_view, tmp_WF_view, tmp_Tst_view, tmp_D, floris, wind.shear)
        if !isnothing(alloc)
            alloc.n += 1
            alloc.floris += a
        end

        a = @allocated begin
            T_red = prod(T_red_arr)
            wf.red_arr[iT, vcat(wf.dep[iT], iT)] = T_red_arr
            T_addedTI = sqrt(sum(T_aTI_arr .^ 2))
            wf.Weight[iT] = T_weight
        end
        if !isnothing(alloc)
            alloc.begin1 += a
        end

        a = @allocated if wf.D[end] <= 0
            # Reuse buffers for distance and plotting calculations
            dists_view = @view ub.dists_buffer[1:(tmp_nT - 1)]
            plot_WF_view = @view ub.plot_WF_buffer[1:(tmp_nT - 1), :]
            plot_OP_view = @view ub.plot_OP_buffer[1:(tmp_nT - 1), :]
            
            dists_view .= 0.0
            plot_WF_view .= 0.0
            plot_OP_view .= 0.0
            
            for iiT in 1:(tmp_nT - 1)
                OP1_i_f, OP1_r, OP2_i_f, OP2_r = wf.intOPs[iT][iiT, :]
                OP1_i = Int(round(OP1_i_f))
                OP2_i = Int(round(OP2_i_f))
                OPi_l = OP1_r * wf.States_OP[OP1_i, :] + OP2_r * wf.States_OP[OP2_i, :]
                plot_OP_view[iiT, :] = OPi_l[1:2]
                plot_WF_view[iiT, :] = OP1_r * wf.States_WF[OP1_i, :] + OP2_r * wf.States_WF[OP2_i, :]
                dists_view[iiT] = norm(OPi_l[1:2] .- wf.posBase[iT,1:2])
            end

            I = sortperm(dists_view)
            if length(I) == 1
                Ufree = plot_WF_view[I[1], 1]
                T_Ueff = T_red * Ufree
            else
                a = plot_OP_view[I[1], :]'
                b = plot_OP_view[I[2], :]'
                c =wf.posBase[iT, 1:2]'
                d = clamp((dot(b - a, c - a)) / dot(b - a, b - a), 0.0, 1.0)
                r1, r2 = 1.0 - d, d
                Ufree = r1 * plot_WF_view[I[1], 1] + r2 * plot_WF_view[I[2], 1]
                T_Ueff = T_red * Ufree
            end
        end
        if !isnothing(alloc)
            alloc.if2 += a
        end

        ub.M_buffer[iT, :] = [T_red, T_addedTI, T_Ueff]

        wS = sum(wf.Weight[iT])
        if wS > 0
           wf.Weight[iT] = wf.Weight[iT] ./ wS
        else
           wf.Weight[iT] .= 0.0
        end
    end
    if !isnothing(alloc)
        alloc.for1 += a
    end

    return ub.M_buffer, wf
end

@with_kw_noshow mutable struct Allocations
    iterateOPs::Int64=0
    perturbationOfTheWF::Int64=0
    findTurbineGroups::Int64=0
    interpolateOPs::Int64=0
    setUpTmpWFAndRun::Int64=0
    correctVel::Int64=0
    correctDir::Int64=0
    correctTI::Int64=0
    getYaw::Int64=0
    getPower::Int64=0
    calcFlowField::Int64=0
    cff_X::Int64=0
    cff_Y::Int64=0
    getMeasurementsP::Int64=0
    gmp_mx::Int64=0
    gmp_buffers::Int64=0
    gmp_alloc2::Int64=0 # allocations of setUpTmpWFAndRun! inside of the threaded loop
end

function Base.show(io::IO, allocs::Allocations)
    println(io, "Allocations:")
    for field_name in fieldnames(typeof(allocs))
        value = getfield(allocs, field_name)
        if value > 5e7
            gb_value = value / 1e9
            println(io, "  $field_name: $(round(gb_value, digits=3)) GB")
        end
    end
end

"""
    runFLORIDyn(plt, set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, vis::Vis,
                floridyn::FloriDyn, floris::Floris; rmt_plot_fn=nothing, msr=VelReduction) -> (WindFarm, DataFrame, Matrix)

Main entry point for the FLORIDyn closed-loop simulation.

# Arguments
- `plt`: Plot object for live visualization during simulation
- `set::Settings`: Simulation settings and configuration parameters.
- `wf::WindFarm`: See: [WindFarm](@ref) simulation state, including turbine and wind farm states.
- `wind::Wind`: See: [Wind](@ref) field settings.
- `sim::Sim`: Simulation state or configuration object. See: [`Sim`](@ref)
- `con::Con`: Controller object or control parameters. See: [`Con`](@ref)
- `vis::Vis`: Visualization settings controlling online plotting and animation. See: [`Vis`](@ref)
- `floridyn::FloriDyn`: Parameters specific to the FLORIDyn model. See: [`FloriDyn`](@ref)
- `floris::Floris`: Parameters specific to the FLORIS model. See: [`Floris`](@ref)

# Keyword Arguments
- `rmt_plot_fn`: Optional remote plotting function for intermediate simulation results. When provided, this function 
  is called remotely (using `@spawnat 2`) to plot flow field visualization on a separate worker process.
  The function should accept parameters `(wf, X, Y, Z, vis, t_rel; msr=VelReduction)` where `wf` is the wind farm state,
  `X`, `Y`, `Z` are flow field coordinates and velocities, `vis` contains visualization settings, and `t_rel` 
  is the relative simulation time. Defaults to `nothing` for local plotting.
- `msr`: Measurement type for velocity reduction calculations. Defaults to `VelReduction`.

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
function runFLORIDyn(plt, set::Settings, wf::WindFarm, wind::Wind, sim::Sim, con::Con, 
                          vis::Vis, floridyn::FloriDyn, floris::Floris; rmt_plot_fn=nothing, 
                          msr=VelReduction, debug=nothing)
    alloc = Allocations()
    nT      = wf.nT
    sim_steps    = sim.n_sim_steps
    ma       = zeros(sim_steps * nT, 6)
    ma[:, 1] .= 1.0  # Set first column to 1
    vm_int   = Vector{Matrix{Float64}}(undef, sim_steps)

    sim_time = sim.start_time
    plot_state = nothing  # Initialize animation state
    
    buffers = FLORIDyn.IterateOPsBuffers(wf)
    # Create unified buffers for all operations
    unified_buffers = create_unified_buffers(wf)
    # Create buffers for interpolateOPs! (will be resized as needed)
    intOPs_buffers = [Matrix{Float64}(undef, 0, 4) for _ in 1:wf.nT]
    
    for it in 1:sim_steps
        sim.sim_step = it

        # ========== PREDICTION ==========
        a = @allocated iterateOPs!(set.iterate_mode, wf, sim, floris, floridyn, buffers)
        alloc.iterateOPs += a

        # ========== Wind Field Perturbation ==========
        a = @allocated perturbationOfTheWF!(wf, wind)
        alloc.perturbationOfTheWF += a

        # ========== Get FLORIS reductions ==========
        a = @allocated wf.dep = findTurbineGroups(wf, floridyn)
        alloc.findTurbineGroups += a
        if sim_steps == 1 && ! isnothing(debug)
            debug[2] = deepcopy(wf)
        end        
        a = @allocated begin
            # Resize buffers if dependencies changed
            for iT in 1:wf.nT
                if size(intOPs_buffers[iT], 1) != length(wf.dep[iT])
                    intOPs_buffers[iT] = zeros(length(wf.dep[iT]), 4)
                end
            end
            wf.intOPs = interpolateOPs!(unified_buffers, intOPs_buffers, wf)
        end
        if sim_steps == 1 && ! isnothing(debug)
            debug[1] = deepcopy(wf)
        end
        alloc.interpolateOPs += a
        a = @allocated tmpM, wf = setUpTmpWFAndRun(set, wf, floris, wind)
        if sim_steps == 1
            # println("intOPs: $(wf.intOPs)")
        end
        alloc.setUpTmpWFAndRun += a

        ma[(it-1)*nT+1 : it*nT, 2:4] .= tmpM
        ma[(it-1)*nT+1 : it*nT, 1]   .= sim_time
        wf.States_T[wf.StartI, 3] = tmpM[:, 2]
   
        vm_int[it] = wf.red_arr

        # ========== wind field corrections ==========
        a = @allocated wf, wind = correctVel(set.cor_vel_mode, set, wf, wind, sim_time, floris, tmpM)
        alloc.correctVel += a
        a = @allocated correctDir!(set.cor_dir_mode, set, wf, wind, sim_time)
        alloc.correctDir += a
        a = @allocated correctTI!(set.cor_turb_mode, set, wf, wind, sim_time)
        alloc.correctTI += a

        # Save free wind speed as measurement
        ma[(it-1)*nT+1 : it*nT, 5] = wf.States_WF[wf.StartI, 1]

        # ========== Get Control settings ==========
        a = @allocated wf.States_T[wf.StartI, 2] = (
            wf.States_WF[wf.StartI, 2] .-
                getYaw(set.control_mode, con.yaw_data, (1:nT), sim_time)'
        )
        alloc.getYaw += a

        # ========== Calculate Power ==========
        P = getPower(wf, tmpM, floris, con)
        ma[(it-1)*nT+1:it*nT, 6] = P

        # ========== Live Plotting ============
        a = @allocated if vis.online
            t_rel = sim_time-sim.start_time
            if mod(t_rel, vis.up_int) == 0
                Z, X, Y = calcFlowField(set, wf, wind, floris; plt, alloc)
                if isnothing(rmt_plot_fn)
                    plot_state = plotFlowField(plot_state, plt, wf, X, Y, Z, vis, t_rel; msr)
                    plt.pause(0.01)
                else
                    # @info "time: $t_rel, plotting with rmt_plot_fn"
                    @spawnat 2 rmt_plot_fn(wf, X, Y, Z, vis, t_rel; msr=msr)
                end
            end
        end
        alloc.calcFlowField += a

        sim_time += sim.time_step
    end
    # Convert `ma` to DataFrame and scale measurements
    md = DataFrame(
        (ma * diagm([1; 100; 100; 1; 1; 1e-6])),
        [:Time, :ForeignReduction, :AddedTurbulence, :EffWindSpeed, :FreeWindSpeed, :PowerGen]
    )
    mi = hcat(md.Time, hcat(vm_int...)')
    println("Memory allocations: $(alloc)")
    return wf, md, mi
end
